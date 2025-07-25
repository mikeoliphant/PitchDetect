using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PitchDetect
{
    public class EnvelopeFollower
    {
        float attackCoef;
        float decayCoef;
        float envelope = 0;

        public EnvelopeFollower(int sampleRate, float attackMS, float decayMS)
        {
            attackCoef = (float)Math.Exp(Math.Log(0.01) / (attackMS * sampleRate * 0.001f));
            decayCoef = (float)Math.Exp(Math.Log(0.01) / (decayMS * sampleRate * 0.001f));
        }

        public float Process(ReadOnlySpan<float> samples)
        {
            for (int i = 0; i < samples.Length; i++)
            {
                float tmp = Math.Abs(samples[i]);

                if (tmp > envelope)
                    envelope = attackCoef * (envelope - tmp) + tmp;
                else
                    envelope = decayCoef * (envelope - tmp) + tmp;
            }

            return envelope;
        }
    }

    public class PeakFollower
    {
        float max = 0;
        float peak = 0;
        int sumSamples;
        int samplesSoFar = 0;

        public PeakFollower(int sampleRate)
        {
            sumSamples = sampleRate / 40;
        }

        public void Reset()
        {
            max = 0;
            peak = 0;
            samplesSoFar = 0;
        }

        public float Process(ReadOnlySpan<float> samples)
        {
            for (int i = 0; i < samples.Length; i++)
            {
                float amp = Math.Abs(samples[i]);

                max = Math.Max(amp, max);

                samplesSoFar++;

                if (samplesSoFar == sumSamples)
                {
                    peak = max;

                    samplesSoFar = 0;
                    max = 0;
                }
            }

            return Math.Max(peak, max);
        }
    }
}
