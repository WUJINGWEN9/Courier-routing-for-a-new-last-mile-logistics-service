using System;
using ILOG.Concert;
using ILOG.CPLEX;
using System.IO;
using System.Numerics;

namespace 帮我买最大化利润
{
    class Program
    {
        const int R = 50, S = 30, K = 60, O = 50, W = 20, MM = 9999, M = 1000, MaxN = 3, NN = 50, SS = 167, DD = 241; 

        const double t0 = 12, f = 2;

        static void Main(string[] args)
        {
            Program a = new Program();
            Random random = new Random();

            #region 参数

            string[] strArray1 = File.ReadAllLines(@"C:\Users\zhenlu\Desktop\实验\试验distance.txt");
            double[][] distance = new double[SS + 1][];
            int ee = 0;
            for (int i = 1; i < SS + 1; i++)
            {
                distance[i] = new double[DD + 1];
                for (int j = 1; j < DD + 1; j++)
                {
                    distance[i][j] = Convert.ToDouble(strArray1[ee].ToString());
                    ee++;
                }
            }
            //for (int i = 1; i < SS + 1; i++)
            //{
            //    for (int j = 1; j < DD + 1; j++)
            //    {
            //        Console.WriteLine("distance[i{0}][j{1}]={2}", i, j, distance[i][j]);
            //    }
            //}
            int snode = 0;
            int rnode = 0;
            double[][] d_ij = new double[S + 1][];
            for (int i = 1; i < S + 1; i++)
            {
                d_ij[i] = new double[R + 1];
                for (int j = 1; j < R + 1; j++)
                {
                    d_ij[i][j] = 0;
                }
            }
            for (int s = 1; s < S + 1; s++)
            {
                for (int r = 1; r < R + 1; r++)
                {
                    snode = random.Next(1, SS + 1);
                    rnode = random.Next(1, DD + 1);
                    d_ij[s][r] = distance[snode][rnode];
                }
            }

            #region p_rs[r][s] 商家s出售已知订单r对应商品的price
            double[][] p_rs = new double[R + 1][];
            double n = 0.2;
            for (int r = 1; r < R + 1; r++)
            {
                p_rs[r] = new double[S + 1];
                for (int s = 1; s < S + 1; s++)
                {
                    p_rs[r][s] = random.Next(5, 20);
                    //Console.WriteLine(" 潜在损失p_rs[r{0}][s{1}]={2}", r, s, p_rs[r][s]);
                }
            }
            #endregion

            #region d_ijk[i][j][k] 代购员k从节点i到节点j之间所需的平均行驶距离 0代表配送员起始位置
            double[][][] d_ijk = new double[S + 1][][];
            for (int i = 0; i < S + 1; i++)
            {
                d_ijk[i] = new double[R + 1][];
                for (int j = 1; j < R + 1; j++)
                {
                    d_ijk[i][j] = new double[K];
                    for (int k = 0; k < K; k++)
                    {
                        if (i > 0)
                        {
                            d_ijk[i][j][k] = Math.Round(d_ij[i][j], 2);
                        }
                    }
                }
            }
            for (int s = 1; s < S + 1; s++)
            {
                for (int k = 0; k < K; k++)
                {
                    d_ijk[0][s][k] = Math.Round(random.Next(3) + random.NextDouble(), 2);
                }
            }

            #endregion

            #region t_ijk[i][j][k] 代购员k从节点i到节点j之间所需的平均行驶时间 速度40km/h
            double[][][] t_ijk = new double[S + 1][][];
            for (int i = 0; i < S + 1; i++)
            {
                t_ijk[i] = new double[R + 1][];
                for (int j = 1; j < R + 1; j++)
                {
                    t_ijk[i][j] = new double[K];
                    for (int k = 0; k < K; k++)
                    {
                        t_ijk[i][j][k] = Math.Round(d_ijk[i][j][k] / 0.67, 2);
                    }
                }
            }

            //for (int i = 0; i < S + 1; i++)
            //{
            //    for (int j = 1; j < R + 1; j++)
            //    {
            //        for (int k = 0; k < K; k++)
            //        {
            //            Console.WriteLine("d_ijk[i{0}][j{1}][k{2}]={3}  t={4}", i, j, k, d_ijk[i][j][k], t_ijk[i][j][k]);
            //        }
            //    }
            //}
            #endregion

            #region d_ijk[i][j][k] 代购员k从节点i到节点j之间所需的平均行驶成本 0代表配送员起始位置
            double m = 0.2;
            for (int i = 0; i < S + 1; i++)
            {
                for (int j = 1; j < R + 1; j++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        d_ijk[i][j][k] = d_ijk[i][j][k] * m;
                        //Console.WriteLine(" d_ijk[i{0}][j{1}][k{2}]={3}", i, j, k, d_ijk[i][j][k]);
                    }
                }
            }

            #endregion

            #region z_rs[r][s] 代购员在提供订单r对应商品的商家s处的等待时间
            double[][] z_rs = new double[R + 1][];
            for (int r = 1; r < R + 1; r++)
            {
                z_rs[r] = new double[S + 1];
                for (int s = 1; s < S + 1; s++)
                {
                    z_rs[r][s] = random.Next(1, 5);
                    //Console.WriteLine("z_rs[r{0}][s{1}]={2}", r, s, z_rs[r][s]);
                }
            }
            #endregion

            #region u_r[r] 已知订单r的最晚可送达时间
            double[] u_r = new double[R + 1];
            for (int r = 1; r < R + 1; r++)
            {
                u_r[r] = t0 + 60;
            }
            #endregion

            #region 二阶段参数

            #region d_wijk[w][i][j][k] 情景ω下代购员k从节点i到节点j之间所需的平均行驶距离 0代表配送员起始位置
            double[][][][] d_wijk = new double[W][][][];
            for (int w = 0; w < W; w++)
            {
                d_wijk[w] = new double[R + 1][][];
                for (int i = 0; i < R + 1; i++)
                {
                    d_wijk[w][i] = new double[O + 1][];
                    for (int j = 1; j < O + 1; j++)
                    {
                        d_wijk[w][i][j] = new double[K];
                        for (int k = 0; k < K; k++)
                        {
                            d_wijk[w][i][j][k] = Math.Round(random.Next(3) + random.NextDouble(), 2) * m;
                        }
                    }
                }
            }

            //for (int w = 0; w < W; w++)
            //{
            //    for (int i = 0; i < R + 1; i++)
            //    {
            //        for (int j = 1; j < O + 1; j++)
            //        {
            //            for (int k = 0; k < K; k++)
            //            {
            //                Console.WriteLine("d_wijk[w{0}][i{1}][j{2}][k{3}]={4}", w, i, j, k, d_wijk[w][i][j][k]);
            //            }
            //        }
            //    }
            //}

            #endregion

            #region pai_w[w] 场景ω出现的概率
            double[] pai_w = new double[W];
            for (int w = 0; w < W; w++)
            {
                pai_w[w] = 0.05;
                Console.WriteLine("pai_w[w{0}]={1}", w, pai_w[w]);
            }
            #endregion

            #endregion

            #endregion

        }
    }
}
