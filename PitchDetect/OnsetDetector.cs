using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PitchDetect
{
    public class EnvelopeFollower
    {
        float attackCoefficient;
        float releaseCoefficient;
        float env = 0;

        public EnvelopeFollower(float attackMS, float releaseMS, int sampleRate)
        {
            attackCoefficient = MathF.Exp(-1.0f / (attackMS * 0.001f * (float)sampleRate));
            releaseCoefficient = MathF.Exp(-1.0f / (releaseMS * 0.001f * (float)sampleRate));
        }

        public void Reset()
        {
            env = 0;
        }

        public void SetValue(float absValue)
        {
            env = absValue;
        }

        public float Advance(float absValue)
        {
            if (absValue > env)
                env = attackCoefficient * (env - absValue) + absValue;
            else
                env = releaseCoefficient * (env - absValue) + absValue;

            return env;
        }
    }

    public class OnsetDetector
    {
        public float AbsoluteThreshold { get; set; } = 0.05f;
        public float RelativeThreshold { get; set; } = 0.4f;
        public int ResetSamples { get; set; } = 4096;
        public float CurrentMax { get { return env; } }

        int resetSamplesSoFar = 0;
        EnvelopeFollower fastFollower;
        EnvelopeFollower slowFollower;
        float env = 0;
        public float envSlow = 0;

        public OnsetDetector(int sampleRate)
        {
            fastFollower = new(0, 35, sampleRate);
            slowFollower = new(5, 100, sampleRate);
        }

        public void Reset()
        {
            fastFollower.Reset();
            slowFollower.Reset();
        }

        public bool? Detect(ReadOnlySpan<float> audioData)
        {
            float abs = 0;

            bool? detect = null;

            for (int pos = 0; pos < audioData.Length ; pos++)
            {
                abs = MathF.Abs(audioData[pos]);

                env = fastFollower.Advance(abs);
                envSlow = slowFollower.Advance(abs);

                if (resetSamplesSoFar > 0)
                {
                    resetSamplesSoFar--;
                }
                else if (!detect.HasValue)
                {
                    float delta = env - envSlow;

                    if (delta > 0)
                    {
                        if ((env > AbsoluteThreshold) && ((delta / env) > RelativeThreshold))
                        {
                            detect = true;

                            resetSamplesSoFar = ResetSamples;
                        }
                    }
                    else
                    {
                        if ((-delta / env) > RelativeThreshold)
                            detect = false;
                    }

                    if (detect.HasValue)
                    {
                        slowFollower.SetValue(env);
                    }
                }
            }

            return detect;
        }
    }

}
