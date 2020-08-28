using System;
using System.IO;
namespace ConsoleDynamics
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
                Console.Write(" type a number as initial slope for shooting \n type 'simu' to simulate 50 steps of dynamics \n type 'export' to export simulation to csv file\n type 'exit' to quit \n ");/*\label{line:user_interface}*/
                string input = Console.ReadLine();
                if (input.Equals("exit"))
                {
                    break;
                }
                else if (input.Equals("simu"))
                {
                    simulate();
                }
                else if (input.Equals("export"))
                {
                    simulation_export_csv(last_t_index);
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








        static double maxx = 50;
        static double dx = 0.2;
        static int xsteps = (int)Math.Round(maxx / dx);

        static double maxt = 500;
        static double dt = 0.1;
        static int tsteps = (int)Math.Round(maxt / dt);
        static int last_t_index = 0;

        static double friction_factor = 1.01;/*\label{line:friction_factor}*/

        static double B1_maxx = 30;
        static int B1_xsteps = (int)Math.Round(B1_maxx / dx);
        static double[,,] B1_fi;
        static double[,,] B1_fi_rotated;

        static double[,,,] fi;
        static double[,,,,] dfi;
        static double[,,,,] ddfi;


        static void simulate()
        {
            if (fi == null)
            {/*\label{line:initial_condition_if}*/
                fi = new double[3, tsteps, xsteps, xsteps];
                dfi = new double[3, 3, tsteps, xsteps, xsteps];
                ddfi = new double[3, 3, 3, xsteps, xsteps];

                B1_fi = build_B1_particle(0);
                B1_fi_rotated = build_B1_particle(Math.PI);

                int p1x = 30;
                int p1y = 0;

                int p2x = 60;
                int p2y = 80;

                int p3x = 0;
                int p3y = 80;

                for (int ti = 0; ti < 2; ti++)
                {
                    for (int xi = 0; xi < xsteps; xi++)
                    {
                        for (int yi = 0; yi < xsteps; yi++)
                        {
                            //boundary condition
                            fi[0, ti, xi, yi] = 0;
                            fi[1, ti, xi, yi] = 0;
                            fi[2, ti, xi, yi] = 1;
                        }
                    }

                    for (int xi = p1x; xi < p1x + B1_xsteps; xi++)
                    {
                        for (int yi = p1y; yi < p1y + B1_xsteps; yi++)
                        {
                            fi[0, ti, xi, yi] += B1_fi[0, xi - p1x, yi - p1y];
                            fi[1, ti, xi, yi] += B1_fi[1, xi - p1x, yi - p1y];
                            fi[2, ti, xi, yi] += B1_fi[2, xi - p1x, yi - p1y] - 1;
                        }
                    }

                    for (int xi = p2x; xi < p2x + B1_xsteps; xi++)
                    {
                        for (int yi = p2y; yi < p2y + B1_xsteps; yi++)
                        {
                            fi[0, ti, xi, yi] += B1_fi[0, xi - p2x, yi - p2y];
                            fi[1, ti, xi, yi] += B1_fi[1, xi - p2x, yi - p2y];
                            fi[2, ti, xi, yi] += B1_fi[2, xi - p2x, yi - p2y] - 1;
                        }
                    }

                    for (int xi = p3x; xi < p3x + B1_xsteps; xi++)
                    {
                        for (int yi = p3y; yi < p3y + B1_xsteps; yi++)
                        {
                            fi[0, ti, xi, yi] += B1_fi_rotated[0, xi - p3x, yi - p3y];
                            fi[1, ti, xi, yi] += B1_fi_rotated[1, xi - p3x, yi - p3y];
                            fi[2, ti, xi, yi] += B1_fi_rotated[2, xi - p3x, yi - p3y] - 1;
                        }
                    }
                }
            }

            for (int ti = 0; ti < 50; ti++)
            {
                int new_t_index = last_t_index + 1;
                if (new_t_index >= tsteps - 1)
                {
                    Console.WriteLine("end of simulation");
                    break;
                }
                try
                {
                    step(new_t_index);
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    break;
                }
                last_t_index = new_t_index;
            }
        }

        static double[,,] build_B1_particle(double deltaTeta)
        {
            double[,,] fi_out = new double[3, B1_xsteps, B1_xsteps];
            for (int xi = 0; xi < B1_xsteps; xi++)
            {
                double x = -B1_maxx / 2 + xi * dx;
                for (int yi = 0; yi < B1_xsteps; yi++)
                {
                    double y = -B1_maxx / 2 + yi * dx;
                    double r = Math.Sqrt(x * x + y * y);
                    int ri = (int)(Math.Round(r / dr));

                    double teta = Math.Atan2(y, x);
                    fi_out[0, xi, yi] = Math.Sin(profile_function[ri]) * Math.Cos(teta - deltaTeta);
                    fi_out[1, xi, yi] = Math.Sin(profile_function[ri]) * Math.Sin(teta - deltaTeta);
                    fi_out[2, xi, yi] = Math.Cos(profile_function[ri]);
                }
            }
            return fi_out;
        }

        static double A;
        static double BB;
        static double C;
        static double D;
        static double E;
        static double F;
        static double G;
        static double H;
        static double M;
        static double T;
        static double U;
        static double V;
        static double W;
        static double minimum_acceptable_accuracy_parameter = 0.2;

        static void step(int ti)
        {
            Console.WriteLine("simulating t = " + ti);
            calculate_dfi(ti);

            for (int xi = 0; xi < xsteps; xi++)
            {
                for (int yi = 0; yi < xsteps; yi++)
                {
                    fi[0, ti + 1, xi, yi] = 0;
                    fi[1, ti + 1, xi, yi] = 0;
                    fi[2, ti + 1, xi, yi] = 0;

                    double alpha = 1;
                    double[] laplace = new double[3];

                    for (int nu = 0; nu < 3; nu++)
                    {
                        for (int fz = 0; fz < 3; fz++)
                        {
                            alpha += -k * k * raise_index(nu) * dfi[nu, fz, ti, xi, yi] * dfi[nu, fz, ti, xi, yi];
                        }
                    }

                    double beta = alpha;
                    for (int fz = 0; fz < 3; fz++)
                    {
                        beta += k * k * dfi[0, fz, ti, xi, yi] * dfi[0, fz, ti, xi, yi];
                    }

                    for (int fz = 0; fz < 3; fz++)
                    {
                        for (int nu = 0; nu < 3; nu++)
                        {
                            ddfi[1, nu, fz, xi, yi] = (-1.0 / 12 * dfi[nu, fz, ti, (xi + 2) % xsteps, yi] + 2.0 / 3 * dfi[nu, fz, ti, (xi + 1) % xsteps, yi] - 2.0 / 3 * dfi[nu, fz, ti, (xi - 1 + xsteps) % xsteps, yi] + 1.0 / 12 * dfi[nu, fz, ti, (xi - 2 + xsteps) % xsteps, yi]) / dx;
                            ddfi[2, nu, fz, xi, yi] = (-1.0 / 12 * dfi[nu, fz, ti, xi, (yi + 2) % xsteps] + 2.0 / 3 * dfi[nu, fz, ti, xi, (yi + 1) % xsteps] - 2.0 / 3 * dfi[nu, fz, ti, xi, (yi - 1 + xsteps) % xsteps] + 1.0 / 12 * dfi[nu, fz, ti, xi, (yi - 2 + xsteps) % xsteps]) / dx;
                        }
                        ddfi[1, 2, fz, xi, yi] = ddfi[2, 1, fz, xi, yi];

                        laplace[fz] = ddfi[1, 1, fz, xi, yi] + ddfi[2, 2, fz, xi, yi];
                    }

                    double[] S = new double[3];

                    for (int q = 0; q < 3; q++)
                    {
                        S[q] = alpha * laplace[q] - d_potential(q);
                        for (int i = 1; i < 3; i++)
                        {
                            for (int fj = 0; fj < 3; fj++)
                            {
                                S[q] += -k * k * dfi[0, q, ti, xi, yi] * dfi[i, fj, ti, xi, yi] * ddfi[i, 0, fj, xi, yi] + k * k * ddfi[i, 0, q, xi, yi] * dfi[0, fj, ti, xi, yi] * dfi[i, fj, ti, xi, yi];
                                for (int nu = 0; nu < 3; nu++)
                                {
                                    S[q] += -k * k * raise_index(nu) * dfi[i, q, ti, xi, yi] * dfi[nu, fj, ti, xi, yi] * ddfi[i, nu, fj, xi, yi] + k * k * raise_index(nu) * ddfi[i, nu, q, xi, yi] * dfi[i, fj, ti, xi, yi] * dfi[nu, fj, ti, xi, yi] + k * k * raise_index(nu) * dfi[nu, q, ti, xi, yi] * dfi[nu, fj, ti, xi, yi] * ddfi[i, i, fj, xi, yi];
                                }
                            }
                        }
                    }

                    A = beta * fi[2, ti, xi, yi];
                    BB = 0;
                    C = -1 * beta * fi[0, ti, xi, yi];
                    D = S[0] * fi[2, ti, xi, yi] - S[2] * fi[0, ti, xi, yi];
                    E = 0;
                    F = beta * fi[2, ti, xi, yi];
                    G = -1 * beta * fi[1, ti, xi, yi];
                    H = S[1] * fi[2, ti, xi, yi] - S[2] * fi[1, ti, xi, yi];
                    for (int i = 1; i < 3; i++)
                    {
                        A += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 0, ti, xi, yi] + fi[0, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);
                        BB += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 1, ti, xi, yi] + fi[0, ti, xi, yi] * dfi[i, 1, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);
                        C += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 2, ti, xi, yi] + fi[0, ti, xi, yi] * dfi[i, 2, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);

                        E += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 1, ti, xi, yi] + fi[1, ti, xi, yi] * dfi[i, 0, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);
                        F += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 1, ti, xi, yi] * dfi[i, 1, ti, xi, yi] + fi[1, ti, xi, yi] * dfi[i, 1, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);
                        G += k * k * (-1 * fi[2, ti, xi, yi] * dfi[i, 2, ti, xi, yi] * dfi[i, 1, ti, xi, yi] + fi[1, ti, xi, yi] * dfi[i, 2, ti, xi, yi] * dfi[i, 2, ti, xi, yi]);
                    }

                    M = D * dt * dt / 2 + A * fi[0, ti, xi, yi] + BB * fi[1, ti, xi, yi] + C * fi[2, ti, xi, yi] + dt * (A * dfi[0, 0, ti, xi, yi] + BB * dfi[0, 1, ti, xi, yi] + C * dfi[0, 2, ti, xi, yi]);
                    T = H * dt * dt / 2 + E * fi[0, ti, xi, yi] + F * fi[1, ti, xi, yi] + G * fi[2, ti, xi, yi] + dt * (E * dfi[0, 0, ti, xi, yi] + F * dfi[0, 1, ti, xi, yi] + G * dfi[0, 2, ti, xi, yi]);

                    bool solved = solve_based_on_f3(ti, xi, yi);
                    if (!solved)
                    {
                        solved = solve_based_on_f1(ti, xi, yi);
                        if (!solved)
                        {
                            solved = solve_based_on_f2(ti, xi, yi);

                            if (!solved)
                            {
                                throw new Exception("all 2nd order equations have accuracy_parameter lower than minimum_acceptable_accuracy_parameter");
                            }
                        }
                    }
                }
            }
        }

        static void calculate_dfi(int ti)
        {
            for (int xi = 0; xi < xsteps; xi++)
            {
                for (int yi = 0; yi < xsteps; yi++)
                {
                    for (int fz = 0; fz < 3; fz++)
                    {
                        dfi[0, fz, ti, xi, yi] = (fi[fz, ti, xi, yi] - fi[fz, ti - 1, xi, yi]) / dt / friction_factor;/*\label{line:friction_factor_use}*/
                        dfi[1, fz, ti, xi, yi] = (-1 * fi[fz, ti, (xi - 3 + xsteps) % xsteps, yi] + 9 * fi[fz, ti, (xi - 2 + xsteps) % xsteps, yi] - 45 * fi[fz, ti, (xi - 1 + xsteps) % xsteps, yi] + 45 * fi[fz, ti, (xi + 1) % xsteps, yi] - 9 * fi[fz, ti, (xi + 2) % xsteps, yi] + 1 * fi[fz, ti, (xi + 3) % xsteps, yi]) / (60 * dx);
                        dfi[2, fz, ti, xi, yi] = (-1 * fi[fz, ti, xi, (yi - 3 + xsteps) % xsteps] + 9 * fi[fz, ti, xi, (yi - 2 + xsteps) % xsteps] - 45 * fi[fz, ti, xi, (yi - 1 + xsteps) % xsteps] + 45 * fi[fz, ti, xi, (yi + 1) % xsteps] - 9 * fi[fz, ti, xi, (yi + 2) % xsteps] + 1 * fi[fz, ti, xi, (yi + 3) % xsteps]) / (60 * dx);
                    }
                }
            }
        }

        static double d_potential(int i)
        {
            if (i == 0)
            {
                return 0;
            }
            else if (i == 1)
            {
                return 0;
            }
            else
            {
                return m2 * -1;
            }
        }

        static int raise_index(int nu)
        {
            if (nu == 0)
            {
                return 1;
            }
            else
            {
                return -1;
            }
        }

        static bool solve_based_on_f1(int ti, int xi, int yi)
        {
            U = Math.Pow(F * C - G * BB, 2) + Math.Pow(E * C - G * A, 2) + Math.Pow(E * BB - F * A, 2);
            V = -2 * (E * C - G * A) * (T * C - G * M) - 2 * (E * BB - F * A) * (T * BB - F * M);
            W = Math.Pow(T * C - G * M, 2) + Math.Pow(T * BB - F * M, 2) - Math.Pow(G * BB - F * C, 2);

            double delta1 = Math.Pow(V, 2) - 4 * U * W;
            if (delta1 < 0 || Math.Abs(Math.Sqrt(delta1) / (2 * U)) < minimum_acceptable_accuracy_parameter)
            {/*\label{line:delta_checking}*/
                return false;
            }
            double fi1_1 = (-1 * V + Math.Sqrt(delta1)) / (2 * U);
            double fi1_2 = (-1 * V - Math.Sqrt(delta1)) / (2 * U);

            if (Math.Abs(fi1_1 - fi[0, ti, xi, yi]) < Math.Abs(fi1_2 - fi[0, ti, xi, yi]))
            {/*\label{line:better_answer_checking}*/
                fi[0, ti + 1, xi, yi] = fi1_1;
            }
            else
            {
                fi[0, ti + 1, xi, yi] = fi1_2;
            }

            fi[1, ti + 1, xi, yi] = (T * C - G * M - (C * E - G * A) * fi[0, ti + 1, xi, yi]) / (C * F - G * BB);
            fi[2, ti + 1, xi, yi] = (T * BB - F * M - (E * BB - F * A) * fi[0, ti + 1, xi, yi]) / (G * BB - F * C);

            return true;
        }

        static bool solve_based_on_f2(int ti, int xi, int yi)
        {
            U = Math.Pow(F * C - G * BB, 2) + Math.Pow(C * E - A * G, 2) + Math.Pow(BB * E - A * F, 2);
            V = -2 * (T * C - G * M) * (F * C - G * BB) - 2 * (M * E - A * T) * (BB * E - A * F);
            W = Math.Pow(T * C - G * M, 2) + Math.Pow(M * E - A * T, 2) - Math.Pow(C * E - A * G, 2);

            double delta2 = Math.Pow(V, 2) - 4 * U * W;
            if (delta2 < 0 || Math.Abs(Math.Sqrt(delta2) / (2 * U)) < minimum_acceptable_accuracy_parameter)
            {
                return false;
            }
            double fi2_1 = (-1 * V + Math.Sqrt(delta2)) / (2 * U);
            double fi2_2 = (-1 * V - Math.Sqrt(delta2)) / (2 * U);

            if (Math.Abs(fi2_1 - fi[1, ti, xi, yi]) < Math.Abs(fi2_2 - fi[1, ti, xi, yi]))
            {
                fi[1, ti + 1, xi, yi] = fi2_1;
            }
            else
            {
                fi[1, ti + 1, xi, yi] = fi2_2;
            }

            fi[0, ti + 1, xi, yi] = (T * C - G * M - (C * F - G * BB) * fi[1, ti + 1, xi, yi]) / (C * E - G * A);
            fi[2, ti + 1, xi, yi] = (T * A - E * M - (F * A - E * BB) * fi[1, ti + 1, xi, yi]) / (G * A - E * C);

            return true;
        }

        static bool solve_based_on_f3(int ti, int xi, int yi)
        {
            U = Math.Pow(G * BB - F * C, 2) + Math.Pow(C * E - A * G, 2) + Math.Pow(E * BB - F * A, 2);
            V = -2 * (T * BB - F * M) * (G * BB - F * C) - 2 * (M * E - A * T) * (C * E - A * G);
            W = Math.Pow(T * BB - F * M, 2) + Math.Pow(M * E - A * T, 2) - Math.Pow(E * BB - F * A, 2);

            double delta3 = Math.Pow(V, 2) - 4 * U * W;
            if (delta3 < 0 || Math.Abs(Math.Sqrt(delta3) / (2 * U)) < minimum_acceptable_accuracy_parameter)
            {
                return false;
            }
            double fi3_1 = (-1 * V + Math.Sqrt(delta3)) / (2 * U);
            double fi3_2 = (-1 * V - Math.Sqrt(delta3)) / (2 * U);

            if (Math.Abs(fi3_1 - fi[2, ti, xi, yi]) < Math.Abs(fi3_2 - fi[2, ti, xi, yi]))
            {
                fi[2, ti + 1, xi, yi] = fi3_1;
            }
            else
            {
                fi[2, ti + 1, xi, yi] = fi3_2;
            }

            fi[0, ti + 1, xi, yi] = (T * BB - F * M - (G * BB - F * C) * fi[2, ti + 1, xi, yi]) / (E * BB - F * A);
            fi[1, ti + 1, xi, yi] = (T * A - E * M - (G * A - E * C) * fi[2, ti + 1, xi, yi]) / (F * A - E * BB);

            return true;
        }

        static void simulation_export_csv(int ti)
        {
            for (int fi_index = 0; fi_index < 3; fi_index++)
            {
                string file_name = "f" + fi_index + ".csv";
                try
                {
                    File.Delete(file_name);
                    using (FileStream fs = File.OpenWrite(file_name))
                    {
                        StreamWriter sw = new StreamWriter(fs);
                        for (int i = 0; i < xsteps; i++)
                        {
                            string line = "";
                            for (int j = 0; j < xsteps; j++)
                            {
                                line += string.Format("{0:0.0000}", fi[fi_index, ti, i, j]);
                                if (j != xsteps - 1)
                                {
                                    line += ",";
                                }
                            }
                            sw.WriteLine(line);
                        }
                        sw.Close();
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine("file " + file_name + " is not writable");
                }
            }
        }


    }
}
