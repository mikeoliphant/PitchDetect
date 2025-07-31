using System.Numerics;

namespace PitchDetect
{
    public class PitchDetector
    {
        public int MinFrequency { get; set; } = 30;
        public float Threshold { get; set; } = 0.5f;
        public uint SampleRate { get; private set; }
        public ReadOnlySpan<Complex> Correlation { get { return complexData; } }

        Complex[] complexData = null;
        FftFlat.FastFourierTransform fft = null;
        uint windowSize;
        (float Freq, float Corr)[] peaks = new (float Freq, float Corr)[16];

        public PitchDetector(uint sampleRate, uint windowSize, bool zeroPad)
        {
            this.SampleRate = sampleRate;
            this.windowSize = windowSize;

            complexData = new Complex[zeroPad ? windowSize * 2 : windowSize];
            fft = new(complexData.Length);
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

            return GetPitch(complexData, MinFrequency, Threshold);
        }

        public static double HannWindow(int n, int frameSize)
        {
            return 0.5 * (1 - Math.Cos((2 * Math.PI * n) / (frameSize - 1)));
        }

        public static float GetSpectralFlux(ReadOnlySpan<float> spectrum1, ReadOnlySpan<float> spectrum2)
        {
            float tot = 0;

            for (int i = 0; i < spectrum1.Length; i++)
            {
                float diff = spectrum1[i] - spectrum2[i];

                tot += diff * diff;
            }

            return tot / spectrum1.Length;
        }

        public void GetSpectrum(ReadOnlySpan<float> audioData, float[] spectrum)
        {
            for (int i = 0; i < audioData.Length; i++)
            {
                complexData[i] = new Complex(audioData[i] * HannWindow(i, (int)windowSize), 0);
            }

            fft.Forward(new Span<Complex>(complexData, 0, audioData.Length));

            for (int i = 0; i < audioData.Length / 2; i++)
            {
                double fft = Math.Abs(complexData[i].Real + complexData[i].Imaginary);
                double fftMirror = Math.Abs(complexData[complexData.Length - i - 1].Real + complexData[complexData.Length - i - 1].Imaginary);

                spectrum[i] = (float)((fft + fftMirror) * (0.5 + (i / (complexData.Length * 2))));
            }
        }

        void AutoCorrelate(Span<Complex> data)
        {
            fft.Forward(complexData);

            for (int i = 0; i < data.Length; i++)
            {
                data[i] *= Complex.Conjugate(data[i]);
            }

            fft.Inverse(data);
        }

        //float Interpolate(float floatBin, params int[] bins)
        //{
        //    double result = 0; // Initialize result

        //    for (int i = 0; i < bins.Length; i++)
        //    {
        //        // Compute individual terms
        //        // of above formula
        //        double term = spectrum[bins[i]];

        //        for (int j = 0; j < bins.Length; j++)
        //        {
        //            if (j != i)
        //                term = term * (floatBin - bins[j]) /
        //                          (double)(bins[i] - bins[j]);
        //        }

        //        // Add current term to result
        //        result += term;
        //    }

        //    return (float)result;
        //}

        //public float GetEnergy(float frequency)
        //{
        //    float bin = ((frequency * (spectrum.Length * 2)) / SampleRate);

        //    int intBin = (int)bin;

        //    return Interpolate(bin, intBin - 1, intBin, intBin + 1);
        //}

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


        //float GetHarmonicEnergy(float baseFreq)
        //{
        //    //float freq = baseFreq * 3;

        //    //int bin = (int)((freq * (spectrum.Length * 2)) / SampleRate);

        //    //return (float)spectrum[bin];

        //    float totEnergy = 0;

        //    for (int harmonic = 0; harmonic < 5; harmonic++)
        //    {
        //        float freq = baseFreq * (1 + harmonic);

        //        totEnergy += GetEnergy(freq);
        //    }

        //    return totEnergy;
        //}


        int GetPeaks(ReadOnlySpan<Complex> corr, int maxBin, float threshold, (float Freq, float Corr)[] results)
        {
            bool haveNegative = false;
            bool isIncreasing = false;
            double lastValue = double.MaxValue;

            int resultPos = 0;

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

                            results[resultPos] = (freq, peakCorr);

                            resultPos++;

                            if (resultPos == results.Length)
                                return resultPos;
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

            return resultPos;
        }

        float GetPitch(ReadOnlySpan<Complex> corr, float minFreq, float threshold)
        {
            int endBin = Math.Min((int)(SampleRate / minFreq), corr.Length);

            Array.Clear(peaks);

            int numPeaks = GetPeaks(corr, endBin, threshold, peaks);

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

            if (numPeaks > 0)
                return peaks.Where(p => p.Corr > (maxPeak * 0.25f)).FirstOrDefault().Freq;

            return 0;
        }
    }
}
