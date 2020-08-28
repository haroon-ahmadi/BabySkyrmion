using System;
namespace ConsoleProfileFunction
{
    class Program
    {
        static double m2 = 0.1;/*\label{line:m2}*/
        static double k = 1;
        static int B = 1;
        static double dr = 0.001;
        static double maxr = 30;
        static int rsteps = (int)Math.Round(maxr / dr);
        static double[] profile_function;
        static double[] d_profile_function;
        static double[] dd_profile_function;
        static double f0 = Math.PI;/*\label{line:f0}*/

        static void Main()
        {
            while (true)
            {
                Console.Write(" type a number as initial slope for shooting \n type 'exit' to quit \n ");
                string input = Console.ReadLine();
                if (input.Equals("exit"))
                {
                    break;
                }
                else
                {
                    double shooting_parameter;
                    if (double.TryParse(input, out shooting_parameter))
                    {
                        shoot(shooting_parameter);
                    }
                    else
                    {
                        Console.Write("Wrong input.");
                    }
                }
            }
        }
        static void shoot(double shooting_parameter)
        {
            profile_function = new double[rsteps];
            d_profile_function = new double[rsteps];
            dd_profile_function = new double[rsteps];
            profile_function[0] = f0;
            d_profile_function[0] = shooting_parameter;

            for (int i = 0; i < rsteps - 1; i++)
            {
                dd_profile_function[i] = G_function(i);
                profile_function[i + 1] = profile_function[i] + dr * d_profile_function[i] + dr * dr * dd_profile_function[i] / 2;
                d_profile_function[i + 1] = d_profile_function[i] + dr * dd_profile_function[i];
                myPrint(i);
            }
        }
        static double G_function(int i)
        {
            if (i == 0)
            {
                return 0;
            }
            double r = i * dr;
            return (
                    m2 * Math.Sin(profile_function[i]) - d_profile_function[i] / r
                    + 0.5 * B * B / r / r * Math.Sin(2 * profile_function[i])
                    - 0.5 * Math.Pow(B * k * d_profile_function[i] / r, 2) * Math.Sin(2 * profile_function[i])
                    + Math.Pow(B * k * Math.Sin(profile_function[i]), 2) * d_profile_function[i] / r / r / r
                ) / (1 + Math.Pow(B * k * Math.Sin(profile_function[i]) / r, 2));
        }
        static void myPrint(int i)
        {
            Console.WriteLine("f({0:0.0000}) = {1}", i * dr, profile_function[i]);
        }
    }
}