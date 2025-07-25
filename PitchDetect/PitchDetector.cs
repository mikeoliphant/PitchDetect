using System.Numerics;

namespace PitchDetect
{
    public class PitchDetector
    {
        public int MinFrequency { get; set; } = 30;
        public float Threshold { get; set; } = 0.5f;
        public int SampleRate { get; private set; }
        public ReadOnlySpan<double> Spectrum { get { return spectrum; } }
        public ReadOnlySpan<Complex> Correlation { get { return complexData; } }

        Complex[] complexData = null;
        double[] spectrum = null;
        FftFlat.FastFourierTransform fft = null;

        public PitchDetector(int sampleRate, int windowSize)
        {
            this.SampleRate = sampleRate;

            complexData = new Complex[windowSize * 2];
            spectrum = new double[windowSize];
            fft = new(windowSize * 2);
        }

        public float GetPitch(ReadOnlySpan<float> audioData)
        {
            for (int i = 0; i < audioData.Length; i++)
            {
                complexData[i] = new Complex(audioData[i], 0);
            }

            for (int i = audioData.Length; i < complexData.Length; i++)
            {
                complexData[i] = 0;
            }

            AutoCorrelate(complexData);

            return GetPitch(complexData, spectrum, MinFrequency, Threshold);
        }

        void AutoCorrelate(Span<Complex> data)
        {
            fft.Forward(data);

            for (int i = 0; i < data.Length / 2; i++)
            {
                double fft = Math.Abs(data[i].Real + data[i].Imaginary);
                double fftMirror = Math.Abs(data[data.Length - i - 1].Real + data[data.Length - i - 1].Imaginary);

                spectrum[i] = (fft + fftMirror) * (0.5 + (i / (data.Length * 2)));
            }

            for (int i = 0; i < data.Length; i++)
            {
                data[i] *= Complex.Conjugate(data[i]);
            }

            fft.Inverse(data);
        }

        float Interpolate(float floatBin, params int[] bins)
        {
            double result = 0; // Initialize result

            for (int i = 0; i < bins.Length; i++)
            {
                // Compute individual terms
                // of above formula
                double term = spectrum[bins[i]];

                for (int j = 0; j < bins.Length; j++)
                {
                    if (j != i)
                        term = term * (floatBin - bins[j]) /
                                  (double)(bins[i] - bins[j]);
                }

                // Add current term to result
                result += term;
            }

            return (float)result;
        }

        public float GetEnergy(float frequency)
        {
            float bin = ((frequency * (spectrum.Length * 2)) / SampleRate);

            int intBin = (int)bin;

            return Interpolate(bin, intBin - 1, intBin, intBin + 1);
        }

        public float GetGoertzelPower(ReadOnlySpan<float> audioData, float freq)
        {
            float s_prev = 0;
            float s_prev2 = 0;
            float coeff, normalizedfreq, power, s;
            int i;

            normalizedfreq = freq / (float)SampleRate;

            coeff = 2 * (float)Math.Cos(2 * Math.PI * normalizedfreq);

            for (i = 0; i < audioData.Length; i++)
            {
                s = audioData[i] + coeff * s_prev - s_prev2;
                s_prev2 = s_prev;
                s_prev = s;
            }

            power = s_prev2 * s_prev2 + s_prev * s_prev - coeff * s_prev * s_prev2;

            return power;
        }


        float GetHarmonicEnergy(float baseFreq)
        {
            //float freq = baseFreq * 3;

            //int bin = (int)((freq * (spectrum.Length * 2)) / SampleRate);

            //return (float)spectrum[bin];

            float totEnergy = 0;

            for (int harmonic = 0; harmonic < 5; harmonic++)
            {
                float freq = baseFreq * (1 + harmonic);

                totEnergy += GetEnergy(freq);
            }

            return totEnergy;
        }


        List<(float Freq, float Corr)> GetPeaks(ReadOnlySpan<Complex> corr, int maxBin, float threshold)
        {
            List<(float Freq, float Corr)> peaks = new();

            bool haveNegative = false;
            bool isIncreasing = false;
            double lastValue = double.MaxValue;

            for (int i = 1; i <= maxBin; i++)
            {
                if (!haveNegative)
                {
                    if (corr[i].Real > 0)
                    {
                        haveNegative = true;
                    }
                    else
                        continue;
                }

                if (isIncreasing)
                {
                    if (corr[i].Real < lastValue)
                    {
                        if (corr[i].Real > threshold)
                        {
                            float freq = 1.0f / ((float)(i - 1) / (float)SampleRate);
                            float peakCorr = (float)corr[i - 1].Real;

                            peaks.Add((freq, peakCorr));
                        }

                        isIncreasing = false;
                    }
                }
                else
                {
                    if ((corr[i].Real > 0) && (corr[i].Real > lastValue))
                    {
                        isIncreasing = true;
                    }
                }

                lastValue = corr[i].Real;
            }

            return peaks;
        }

        float GetPitch(ReadOnlySpan<Complex> corr, ReadOnlySpan<double> spectrum, float minFreq, float threshold)
        {
            int endBin = Math.Min((int)(SampleRate / minFreq), corr.Length);

            var peaks = GetPeaks(corr, endBin, threshold);

            if (peaks.Count == 0)
            {
                return 0;
            }

            float maxPeak = peaks.Select(p => p.Corr).Max();

            //var peakEnergy = new float[peaks.Count];

            //for (int i = 0; i < peaks.Count; i++)
            //{
            //    float freq = 1.0f / (peaks[i] / (float)SampleRate);

            //    peakEnergy[i] = GetHarmonicEnergy(freq);
            //}


            //int max = 0;
            //float maxEnergy = 0;

            //float peakThreash = (float)maxPeak * 0.25f;

            //for (int i = 0; i < peakEnergy.Length; i++)
            //{
            //    if ((corr[peaks[i]].Real > peakThreash) && (peakEnergy[i] > maxEnergy))
            //    {
            //        return (float)SampleRate / (float)peaks[i];

            //        max = i;
            //        maxEnergy = peakEnergy[i];
            //    }
            //}

            //if (maxEnergy > 0)
            //{
            //    return (float)SampleRate / (float)peaks[max];
            //}

            if (peaks.Any())
                return peaks.Where(p => p.Corr > (maxPeak * 0.25f)).FirstOrDefault().Freq;

            return 0;
        }
    }
}
