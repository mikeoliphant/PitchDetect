using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PitchDetect
{
    public class OnsetDetector
    {
        public float AbsoluteThreshold { get; set; } = 0.01f;
        public float RelativeThreshold { get; set; } = 0.5f;
        public int ResetSamples { get; set; } = 4096;
        public float CurrentMax { get { return env; } }

        float currentRunningPeak = 0;
        int resetSamplesSoFar = 0;
        float env;
        float releaseCoefficient = MathF.Exp(-1.0f / (35 * 0.001f * 48000f));
        float lastEnv;

        public OnsetDetector()
        {

        }

        public void Reset()
        {
            env = lastEnv = 0;
            currentRunningPeak = 0;
        }

        public bool Detect(ReadOnlySpan<float> audioData)
        {
            for (int pos = 0; pos < audioData.Length ; pos++)
            {
                float abs = MathF.Abs(audioData[pos]);

                if (abs > env)
                    env = abs;
                else
                    env = releaseCoefficient * (env - abs) + abs;
            }

            bool detect = false;

            if (resetSamplesSoFar > 0)
            {
                resetSamplesSoFar -= audioData.Length;
                currentRunningPeak = env;
            }
            else
            {
                if (env > lastEnv)
                {
                    float delta = (env - currentRunningPeak);

                    if (delta > AbsoluteThreshold)
                    {
                        detect = (delta / env) > RelativeThreshold;

                        if (detect)
                        {
                            currentRunningPeak = env;

                            resetSamplesSoFar = ResetSamples;
                        }
                    }
                }
                else
                {
                    currentRunningPeak = env;
                }
            }

            lastEnv = env;

            return detect;
        }
    }


}
