using System;
using ILOG.Concert;
using ILOG.CPLEX;
using System.IO;
using System.Numerics;

namespace 随机场景
{
    class Program
    {
        const int R = 8, S = 4, K = 15, O = 8, W = 20, MM = 9999, SS = 167, DD = 241; 

        const double t0 = 12, f = 6;

        public static double[,] firstbelta_rk = new double[R + 1, K];
        public void method1(double[][] p_rs, double[][][] d_ijk, double[][][] t_ijk, double[][] z_rs, double[] u_r, double[][][][] d_wijk, double[] pai_w)
        {
            Cplex model = new Cplex();

            #region 决策变量
            INumVar[][] alpha_rs = new INumVar[R + 1][];
            for (int r = 0; r < R + 1; r++)
            {
                alpha_rs[r] = new INumVar[S + 1];
                alpha_rs[r] = model.NumVarArray(S + 1, 0, 1, NumVarType.Bool);
            }

            INumVar[][] belta_rk = new INumVar[R + 1][];
            for (int r = 0; r < R + 1; r++)
            {
                belta_rk[r] = new INumVar[K];
                belta_rk[r] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
            }

            INumVar[][][] theta_rsk = new INumVar[R + 1][][];
            for (int r = 0; r < R + 1; r++)
            {
                theta_rsk[r] = new INumVar[S + 1][];
                for (int s = 0; s < S + 1; s++)
                {
                    theta_rsk[r][s] = new INumVar[K];
                    theta_rsk[r][s] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }
            }

            #region 二阶段决策变量
            INumVar[][][] gama_wok = new INumVar[W][][];
            for (int w = 0; w < W; w++)
            {
                gama_wok[w] = new INumVar[O + 1][];
                for (int o = 1; o < O + 1; o++)
                {
                    gama_wok[w][o] = new INumVar[K];
                    gama_wok[w][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }
            }

            INumVar[][][][] epsilon_wrok = new INumVar[W][][][];
            for (int w = 0; w < W; w++)
            {
                epsilon_wrok[w] = new INumVar[R + 1][][];
                for (int r = 1; r < R + 1; r++)
                {
                    epsilon_wrok[w][r] = new INumVar[O + 1][];
                    for (int o = 1; o < O + 1; o++)
                    {
                        epsilon_wrok[w][r][o] = new INumVar[K];
                        epsilon_wrok[w][r][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                    }
                }
            }

            //INumVar[][][] tau_wok = new INumVar[W][][];
            //for (int w = 0; w < W; w++)
            //{
            //    tau_wok[w] = new INumVar[O + 1][];
            //    for (int o = 1; o < O + 1; o++)
            //    {
            //        tau_wok[w][o] = new INumVar[K];
            //        tau_wok[w][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
            //    }
            //}
            #endregion

            #endregion

            #region 约束
            INumExpr[] expr3 = new INumExpr[20];//约束2
            for (int r = 1; r < R + 1; r++)
            {
                expr3[0] = alpha_rs[0][0];
                expr3[0] = model.Sum(expr3[0], model.Prod(-1, alpha_rs[0][0]));
                for (int k = 0; k < K; k++)
                {
                    expr3[0] = model.Sum(expr3[0], belta_rk[r][k]);
                }
                model.AddLe(expr3[0], 1);
            }

            INumExpr[] expr33 = new INumExpr[20];//约束3
            INumExpr[] expr333 = new INumExpr[20];//约束3
            for (int r = 1; r < R + 1; r++)
            {
                expr33[0] = alpha_rs[0][0];
                expr33[0] = model.Sum(expr33[0], model.Prod(-1, alpha_rs[0][0]));
                expr333[0] = alpha_rs[0][0];
                expr333[0] = model.Sum(expr333[0], model.Prod(-1, alpha_rs[0][0]));
                for (int s = 1; s < S + 1; s++)
                {
                    expr33[0] = model.Sum(expr33[0], alpha_rs[r][s]);
                }
                for (int k = 0; k < K; k++)
                {
                    expr333[0] = model.Sum(expr333[0], belta_rk[r][k]);
                }
                model.AddEq(expr33[0], expr333[0]);
            }



            INumExpr[] expr4 = new INumExpr[20];//约束4
            for (int k = 0; k < K; k++)
            {
                expr4[0] = alpha_rs[0][0];
                expr4[0] = model.Sum(expr4[0], model.Prod(-1, alpha_rs[0][0]));
                for (int r = 1; r < R + 1; r++)
                {
                    expr4[0] = model.Sum(expr4[0], belta_rk[r][k]);
                }
                model.AddLe(expr4[0], 1);
            }

            for (int r = 1; r < R + 1; r++)  //约束5
            {
                for (int s = 1; s < S + 1; s++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        model.AddEq(theta_rsk[r][s][k], model.Max(model.Sum(alpha_rs[r][s], model.Sum(belta_rk[r][k], -1)), 0));
                    }
                }
            }

            INumExpr[] expr6 = new INumExpr[20]; //约束6
            for (int r = 1; r < R + 1; r++)
            {
                for (int k = 0; k < K; k++)
                {
                    expr6[0] = alpha_rs[0][0];
                    expr6[0] = model.Sum(expr6[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int s = 1; s < S + 1; s++)
                    {
                        expr6[0] = model.Sum(expr6[0], model.Prod(theta_rsk[r][s][k], t0 + t_ijk[0][s][k] + z_rs[r][s] + t_ijk[s][r][k]));
                    }
                    model.AddLe(expr6[0], u_r[r]);
                }
            }

            INumExpr[] expr10 = new INumExpr[20];//约束11
            for (int o = 1; o < O + 1; o++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr10[0] = alpha_rs[0][0];
                    expr10[0] = model.Sum(expr10[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int k = 0; k < K; k++)
                    {
                        expr10[0] = model.Sum(expr10[0], gama_wok[w][o][k]);
                    }
                    model.AddLe(expr10[0], 1);
                }
            }

            INumExpr[] expr12 = new INumExpr[20];//13
            INumExpr[] expr122 = new INumExpr[20];
            for (int k = 0; k < K; k++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr12[0] = alpha_rs[0][0];
                    expr12[0] = model.Sum(expr12[0], model.Prod(-1, alpha_rs[0][0]));
                    expr122[0] = alpha_rs[0][0];
                    expr122[0] = model.Sum(expr122[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int o = 1; o < O + 1; o++)
                    {
                        expr12[0] = model.Sum(expr12[0], gama_wok[w][o][k]);
                    }
                    for (int r = 1; r < R + 1; r++)
                    {
                        expr122[0] = model.Sum(expr122[0], belta_rk[r][k]);
                    }
                    model.AddLe(expr12[0], expr122[0]);
                }
            }

            for (int r = 1; r < R + 1; r++)  //约束12
            {
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        for (int w = 0; w < W; w++)
                        {
                            model.AddEq(epsilon_wrok[w][r][o][k], model.Max(model.Sum(belta_rk[r][k], model.Sum(gama_wok[w][o][k], -1)), 0));
                        }
                    }
                }
            }

            #endregion

            #region 目标函数
            INumExpr obj = alpha_rs[0][0];
            obj = model.Sum(obj, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj1 = alpha_rs[0][0];
            obj1 = model.Sum(obj1, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj2 = alpha_rs[0][0];
            obj2 = model.Sum(obj2, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj5 = alpha_rs[0][0];
            obj5 = model.Sum(obj5, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj4 = alpha_rs[0][0];
            obj4 = model.Sum(obj4, model.Prod(-1, alpha_rs[0][0]));
            for (int r = 1; r < R + 1; r++)
            {
                for (int s = 1; s < S + 1; s++)
                {
                    obj = model.Sum(obj, model.Prod(alpha_rs[r][s], p_rs[r][s]));
                }
            }
            for (int r = 1; r < R + 1; r++)
            {
                for (int s = 1; s < S + 1; s++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        obj1 = model.Sum(obj1, model.Prod(theta_rsk[r][s][k], (d_ijk[0][s][k] + d_ijk[s][r][k])));
                    }
                }
            }
            for (int w = 0; w < W; w++)
            {
                for (int r = 1; r < R + 1; r++)
                {
                    for (int o = 1; o < O + 1; o++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            obj2 = model.Sum(obj2, model.Prod(pai_w[w] * d_wijk[w][r][o][k], epsilon_wrok[w][r][o][k]));
                        }
                    }
                }
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        obj5 = model.Sum(obj5, model.Prod(pai_w[w] * f, gama_wok[w][o][k]));
                    }
                }
            }
            obj = model.Sum(obj, model.Prod(-1, obj1));
            obj4 = obj;
            obj2 = model.Sum(obj5, model.Prod(-1, obj2));
            obj = model.Sum(obj, obj2);
            model.AddMaximize(obj);
            #endregion

            //model.SetOut(null);

            #region 输出
            if (model.Solve())
            {
                Console.WriteLine("method1 = {0}   一阶段利润={1},未来订单利润= {2}", Math.Round(model.GetValue(obj), 2), Math.Round(model.GetValue(obj4), 2), Math.Round(model.GetValue(obj2), 2));
            }
            else
            {
                Console.WriteLine("无解");
            }
            #endregion
        }
        public void method2average(ref double method2, double[][] p_rs, double[][][] d_ijk, double[][][] t_ijk, double[][] z_rs, double[] u_r, double[][][] ar_ijk)
        {
            Cplex model = new Cplex();

            #region 决策变量
            INumVar[][] alpha_rs = new INumVar[R + 1][];
            for (int r = 0; r < R + 1; r++)
            {
                alpha_rs[r] = new INumVar[S + 1];
                alpha_rs[r] = model.NumVarArray(S + 1, 0, 1, NumVarType.Bool);
            }

            INumVar[][] belta_rk = new INumVar[R + 1][];
            for (int r = 0; r < R + 1; r++)
            {
                belta_rk[r] = new INumVar[K];
                belta_rk[r] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
            }

            INumVar[][][] theta_rsk = new INumVar[R + 1][][];
            for (int r = 0; r < R + 1; r++)
            {
                theta_rsk[r] = new INumVar[S + 1][];
                for (int s = 0; s < S + 1; s++)
                {
                    theta_rsk[r][s] = new INumVar[K];
                    theta_rsk[r][s] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }
            }

            #region 二阶段决策变量
            INumVar[][] gama_ok = new INumVar[O + 1][];
            for (int o = 1; o < O + 1; o++)
            {
                gama_ok[o] = new INumVar[K];
                gama_ok[o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
            }

            INumVar[][][] epsilon_rok = new INumVar[R + 1][][];
            for (int r = 1; r < R + 1; r++)
            {
                epsilon_rok[r] = new INumVar[O + 1][];
                for (int o = 1; o < O + 1; o++)
                {
                    epsilon_rok[r][o] = new INumVar[K];
                    epsilon_rok[r][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }
            }
            #endregion

            #endregion

            #region 约束
            INumExpr[] expr3 = new INumExpr[20];//约束2
            for (int r = 1; r < R + 1; r++)
            {
                expr3[0] = alpha_rs[0][0];
                expr3[0] = model.Sum(expr3[0], model.Prod(-1, alpha_rs[0][0]));
                for (int k = 0; k < K; k++)
                {
                    expr3[0] = model.Sum(expr3[0], belta_rk[r][k]);
                }
                model.AddLe(expr3[0], 1);
            }

            INumExpr[] expr33 = new INumExpr[20];//约束3
            INumExpr[] expr333 = new INumExpr[20];//约束3
            for (int r = 1; r < R + 1; r++)
            {
                expr33[0] = alpha_rs[0][0];
                expr33[0] = model.Sum(expr33[0], model.Prod(-1, alpha_rs[0][0]));
                expr333[0] = alpha_rs[0][0];
                expr333[0] = model.Sum(expr333[0], model.Prod(-1, alpha_rs[0][0]));
                for (int s = 1; s < S + 1; s++)
                {
                    expr33[0] = model.Sum(expr33[0], alpha_rs[r][s]);
                }
                for (int k = 0; k < K; k++)
                {
                    expr333[0] = model.Sum(expr333[0], belta_rk[r][k]);
                }
                model.AddEq(expr33[0], expr333[0]);
            }



            INumExpr[] expr4 = new INumExpr[20];//约束4
            for (int k = 0; k < K; k++)
            {
                expr4[0] = alpha_rs[0][0];
                expr4[0] = model.Sum(expr4[0], model.Prod(-1, alpha_rs[0][0]));
                for (int r = 1; r < R + 1; r++)
                {
                    expr4[0] = model.Sum(expr4[0], belta_rk[r][k]);
                }
                model.AddLe(expr4[0], 1);
            }

            for (int r = 1; r < R + 1; r++)  //约束5
            {
                for (int s = 1; s < S + 1; s++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        model.AddEq(theta_rsk[r][s][k], model.Max(model.Sum(alpha_rs[r][s], model.Sum(belta_rk[r][k], -1)), 0));
                    }
                }
            }

            INumExpr[] expr6 = new INumExpr[20]; //约束6
            for (int r = 1; r < R + 1; r++)
            {
                for (int k = 0; k < K; k++)
                {
                    expr6[0] = alpha_rs[0][0];
                    expr6[0] = model.Sum(expr6[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int s = 1; s < S + 1; s++)
                    {
                        expr6[0] = model.Sum(expr6[0], model.Prod(theta_rsk[r][s][k], t0 + t_ijk[0][s][k] + z_rs[r][s] + t_ijk[s][r][k]));
                    }
                    model.AddLe(expr6[0], u_r[r]);
                }
            }

            INumExpr[] expr10 = new INumExpr[20];//约束11
            for (int o = 1; o < O + 1; o++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr10[0] = alpha_rs[0][0];
                    expr10[0] = model.Sum(expr10[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int k = 0; k < K; k++)
                    {
                        expr10[0] = model.Sum(expr10[0], gama_ok[o][k]);
                    }
                    model.AddLe(expr10[0], 1);
                }
            }

            INumExpr[] expr12 = new INumExpr[20];//13
            INumExpr[] expr122 = new INumExpr[20];
            for (int k = 0; k < K; k++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr12[0] = alpha_rs[0][0];
                    expr12[0] = model.Sum(expr12[0], model.Prod(-1, alpha_rs[0][0]));
                    expr122[0] = alpha_rs[0][0];
                    expr122[0] = model.Sum(expr122[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int o = 1; o < O + 1; o++)
                    {
                        expr12[0] = model.Sum(expr12[0], gama_ok[o][k]);
                    }
                    for (int r = 1; r < R + 1; r++)
                    {
                        expr122[0] = model.Sum(expr122[0], belta_rk[r][k]);
                    }
                    model.AddLe(expr12[0], expr122[0]);
                }
            }

            for (int r = 1; r < R + 1; r++)  //约束12
            {
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        for (int w = 0; w < W; w++)
                        {
                            model.AddEq(epsilon_rok[r][o][k], model.Max(model.Sum(belta_rk[r][k], model.Sum(gama_ok[o][k], -1)), 0));
                        }
                    }
                }
            }

            #endregion

            #region 目标函数
            INumExpr obj = alpha_rs[0][0];
            obj = model.Sum(obj, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj1 = alpha_rs[0][0];
            obj1 = model.Sum(obj1, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj2 = alpha_rs[0][0];
            obj2 = model.Sum(obj2, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj5 = alpha_rs[0][0];
            obj5 = model.Sum(obj5, model.Prod(-1, alpha_rs[0][0]));
            INumExpr obj4 = alpha_rs[0][0];
            obj4 = model.Sum(obj4, model.Prod(-1, alpha_rs[0][0]));
            for (int r = 1; r < R + 1; r++)
            {
                for (int s = 1; s < S + 1; s++)
                {
                    obj = model.Sum(obj, model.Prod(alpha_rs[r][s], p_rs[r][s]));
                }
            }
            for (int r = 1; r < R + 1; r++)
            {
                for (int s = 1; s < S + 1; s++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        obj1 = model.Sum(obj1, model.Prod(theta_rsk[r][s][k], (d_ijk[0][s][k] + d_ijk[s][r][k])));
                    }
                }
            }
            for (int w = 0; w < W; w++)
            {
                for (int r = 1; r < R + 1; r++)
                {
                    for (int o = 1; o < O + 1; o++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            obj2 = model.Sum(obj2, model.Prod(ar_ijk[r][o][k], epsilon_rok[r][o][k]));
                        }
                    }
                }
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        obj5 = model.Sum(obj5, model.Prod(f, gama_ok[o][k]));
                    }
                }
            }
            obj = model.Sum(obj, model.Prod(-1, obj1));
            obj4 = obj;
            obj2 = model.Sum(obj5, model.Prod(-1, obj2));
            obj = model.Sum(obj, obj2);
            model.AddMaximize(obj);
            #endregion
            model.SetOut(null);
            #region 输出
            if (model.Solve())
            {
                Console.WriteLine("************************ method2一阶段目标值={0}", Math.Round(model.GetValue(obj4), 2));
                method2 = Math.Round(model.GetValue(obj4), 2);
                for (int r = 1; r < R + 1; r++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        firstbelta_rk[r, k] = model.GetValue(belta_rk[r][k]);
                        //Console.WriteLine("firstbelta_rk[r{0}][k{1}]={2}", r, k, firstbelta_rk[r, k]);
                    }
                }
            }
            else
            {
                Console.WriteLine("无解");
            }
            #endregion
        }
        public void method2second(ref double method2, double[][][][] d_wijk, double[] pai_w)
        {
            Cplex model = new Cplex();


            #region 二阶段决策变量
            INumVar[][][] gama_wok = new INumVar[W][][];
            for (int w = 0; w < W; w++)
            {
                gama_wok[w] = new INumVar[O + 1][];
                for (int o = 0; o < O + 1; o++)
                {
                    gama_wok[w][o] = new INumVar[K];
                    gama_wok[w][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }
            }

            INumVar[][][][] epsilon_wrok = new INumVar[W][][][];
            for (int w = 0; w < W; w++)
            {
                epsilon_wrok[w] = new INumVar[R + 1][][];
                for (int r = 1; r < R + 1; r++)
                {
                    epsilon_wrok[w][r] = new INumVar[O + 1][];
                    for (int o = 1; o < O + 1; o++)
                    {
                        epsilon_wrok[w][r][o] = new INumVar[K];
                        epsilon_wrok[w][r][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                    }
                }
            }

            #endregion


            #region 约束           
            INumExpr[] expr10 = new INumExpr[20];//约束10
            for (int o = 1; o < O + 1; o++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr10[0] = gama_wok[0][0][0];
                    expr10[0] = model.Sum(expr10[0], model.Prod(-1, gama_wok[0][0][0]));
                    for (int k = 0; k < K; k++)
                    {
                        expr10[0] = model.Sum(expr10[0], gama_wok[w][o][k]);
                    }
                    model.AddLe(expr10[0], 1);
                }
            }

            INumExpr[] expr12 = new INumExpr[20];//13
            INumExpr[] expr122 = new INumExpr[20];
            for (int k = 0; k < K; k++)
            {
                for (int w = 0; w < W; w++)
                {
                    expr12[0] = gama_wok[0][0][0];
                    expr12[0] = model.Sum(expr12[0], model.Prod(-1, gama_wok[0][0][0]));
                    expr122[0] = gama_wok[0][0][0]; ;
                    expr122[0] = model.Sum(expr122[0], model.Prod(-1, gama_wok[0][0][0]));
                    for (int o = 1; o < O + 1; o++)
                    {
                        expr12[0] = model.Sum(expr12[0], gama_wok[w][o][k]);
                    }
                    for (int r = 1; r < R + 1; r++)
                    {
                        expr122[0] = model.Sum(expr122[0], firstbelta_rk[r, k]);
                    }
                    model.AddLe(expr12[0], expr122[0]);
                }
            }


            for (int r = 1; r < R + 1; r++)  //约束12
            {
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        for (int w = 0; w < W; w++)
                        {
                            model.AddEq(epsilon_wrok[w][r][o][k], model.Max(model.Sum(firstbelta_rk[r, k], model.Sum(gama_wok[w][o][k], -1)), 0));
                        }
                    }
                }
            }
            #endregion

            #region 目标函数
            INumExpr obj = gama_wok[0][0][0];
            obj = model.Sum(obj, model.Prod(-1, gama_wok[0][0][0]));
            INumExpr obj1 = gama_wok[0][0][0];
            obj1 = model.Sum(obj1, model.Prod(-1, gama_wok[0][0][0]));
            for (int w = 0; w < W; w++)
            {
                for (int r = 1; r < R + 1; r++)
                {
                    for (int o = 1; o < O + 1; o++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            obj1 = model.Sum(obj1, model.Prod(pai_w[w] * d_wijk[w][r][o][k], epsilon_wrok[w][r][o][k]));
                        }
                    }
                }
                for (int o = 1; o < O + 1; o++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        obj = model.Sum(obj, model.Prod(pai_w[w] * f, gama_wok[w][o][k]));
                    }
                }
            }
            obj = model.Sum(obj, model.Prod(-1, obj1));
            model.AddMaximize(obj);
            #endregion

            model.SetOut(null);
            #region 输出
            if (model.Solve())
            {
                //Console.WriteLine("二阶段目标值= {0}", Math.Round(model.GetValue(obj), 2));
                method2 = method2 + Math.Round(model.GetValue(obj), 2);
            }
            else
            {
                Console.WriteLine("无解");
            }
            #endregion
        }
        public void method3(ref double method3, double[][] p_rs, double[][][] d_ijk, double[][][] t_ijk, double[][] z_rs, double[] u_r, double[][][][] d_wijk, double[] pai_w)
        {
            double[] sum = new double[W];
            for (int w = 0; w < W; w++)
            {
                Cplex model = new Cplex();

                #region 决策变量
                INumVar[][] alpha_rs = new INumVar[R + 1][];
                for (int r = 0; r < R + 1; r++)
                {
                    alpha_rs[r] = new INumVar[S + 1];
                    alpha_rs[r] = model.NumVarArray(S + 1, 0, 1, NumVarType.Bool);
                }

                INumVar[][] belta_rk = new INumVar[R + 1][];
                for (int r = 0; r < R + 1; r++)
                {
                    belta_rk[r] = new INumVar[K];
                    belta_rk[r] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                }

                INumVar[][][] theta_rsk = new INumVar[R + 1][][];
                for (int r = 0; r < R + 1; r++)
                {
                    theta_rsk[r] = new INumVar[S + 1][];
                    for (int s = 0; s < S + 1; s++)
                    {
                        theta_rsk[r][s] = new INumVar[K];
                        theta_rsk[r][s] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                    }
                }

                #region 二阶段决策变量
                INumVar[][][] gama_wok = new INumVar[W][][];
                //for (int w = 0; w < W; w++)
                {
                    gama_wok[w] = new INumVar[O + 1][];
                    for (int o = 1; o < O + 1; o++)
                    {
                        gama_wok[w][o] = new INumVar[K];
                        gama_wok[w][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                    }
                }

                INumVar[][][][] epsilon_wrok = new INumVar[W][][][];
                //for (int w = 0; w < W; w++)
                {
                    epsilon_wrok[w] = new INumVar[R + 1][][];
                    for (int r = 1; r < R + 1; r++)
                    {
                        epsilon_wrok[w][r] = new INumVar[O + 1][];
                        for (int o = 1; o < O + 1; o++)
                        {
                            epsilon_wrok[w][r][o] = new INumVar[K];
                            epsilon_wrok[w][r][o] = model.NumVarArray(K, 0, 1, NumVarType.Bool);
                        }
                    }
                }

                #endregion

                #endregion

                #region 约束
                INumExpr[] expr3 = new INumExpr[20];//约束2
                for (int r = 1; r < R + 1; r++)
                {
                    expr3[0] = alpha_rs[0][0];
                    expr3[0] = model.Sum(expr3[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int k = 0; k < K; k++)
                    {
                        expr3[0] = model.Sum(expr3[0], belta_rk[r][k]);
                    }
                    model.AddLe(expr3[0], 1);
                }

                INumExpr[] expr33 = new INumExpr[20];//约束3
                INumExpr[] expr333 = new INumExpr[20];//约束3
                for (int r = 1; r < R + 1; r++)
                {
                    expr33[0] = alpha_rs[0][0];
                    expr33[0] = model.Sum(expr33[0], model.Prod(-1, alpha_rs[0][0]));
                    expr333[0] = alpha_rs[0][0];
                    expr333[0] = model.Sum(expr333[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int s = 1; s < S + 1; s++)
                    {
                        expr33[0] = model.Sum(expr33[0], alpha_rs[r][s]);
                    }
                    for (int k = 0; k < K; k++)
                    {
                        expr333[0] = model.Sum(expr333[0], belta_rk[r][k]);
                    }
                    model.AddEq(expr33[0], expr333[0]);
                }



                INumExpr[] expr4 = new INumExpr[20];//约束4
                for (int k = 0; k < K; k++)
                {
                    expr4[0] = alpha_rs[0][0];
                    expr4[0] = model.Sum(expr4[0], model.Prod(-1, alpha_rs[0][0]));
                    for (int r = 1; r < R + 1; r++)
                    {
                        expr4[0] = model.Sum(expr4[0], belta_rk[r][k]);
                    }
                    model.AddLe(expr4[0], 1);
                }

                for (int r = 1; r < R + 1; r++)  //约束5
                {
                    for (int s = 1; s < S + 1; s++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            model.AddEq(theta_rsk[r][s][k], model.Max(model.Sum(alpha_rs[r][s], model.Sum(belta_rk[r][k], -1)), 0));
                        }
                    }
                }

                INumExpr[] expr6 = new INumExpr[20]; //约束6
                for (int r = 1; r < R + 1; r++)
                {
                    for (int k = 0; k < K; k++)
                    {
                        expr6[0] = alpha_rs[0][0];
                        expr6[0] = model.Sum(expr6[0], model.Prod(-1, alpha_rs[0][0]));
                        for (int s = 1; s < S + 1; s++)
                        {
                            expr6[0] = model.Sum(expr6[0], model.Prod(theta_rsk[r][s][k], t0 + t_ijk[0][s][k] + z_rs[r][s] + t_ijk[s][r][k]));
                        }
                        model.AddLe(expr6[0], u_r[r]);
                    }
                }

                INumExpr[] expr10 = new INumExpr[20];//约束11
                for (int o = 1; o < O + 1; o++)
                {
                    //for (int w = 0; w < W; w++)
                    {
                        expr10[0] = alpha_rs[0][0];
                        expr10[0] = model.Sum(expr10[0], model.Prod(-1, alpha_rs[0][0]));
                        for (int k = 0; k < K; k++)
                        {
                            expr10[0] = model.Sum(expr10[0], gama_wok[w][o][k]);
                        }
                        model.AddLe(expr10[0], 1);
                    }
                }

                INumExpr[] expr12 = new INumExpr[20];//13
                INumExpr[] expr122 = new INumExpr[20];
                for (int k = 0; k < K; k++)
                {
                    //for (int w = 0; w < W; w++)
                    {
                        expr12[0] = alpha_rs[0][0];
                        expr12[0] = model.Sum(expr12[0], model.Prod(-1, alpha_rs[0][0]));
                        expr122[0] = alpha_rs[0][0];
                        expr122[0] = model.Sum(expr122[0], model.Prod(-1, alpha_rs[0][0]));
                        for (int o = 1; o < O + 1; o++)
                        {
                            expr12[0] = model.Sum(expr12[0], gama_wok[w][o][k]);
                        }
                        for (int r = 1; r < R + 1; r++)
                        {
                            expr122[0] = model.Sum(expr122[0], belta_rk[r][k]);
                        }
                        model.AddLe(expr12[0], expr122[0]);
                    }
                }

                for (int r = 1; r < R + 1; r++)  //约束12
                {
                    for (int o = 1; o < O + 1; o++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            //for (int w = 0; w < W; w++)
                            {
                                model.AddEq(epsilon_wrok[w][r][o][k], model.Max(model.Sum(belta_rk[r][k], model.Sum(gama_wok[w][o][k], -1)), 0));
                            }
                        }
                    }
                }

                #endregion

                #region 目标函数
                INumExpr obj = alpha_rs[0][0];
                obj = model.Sum(obj, model.Prod(-1, alpha_rs[0][0]));
                INumExpr obj1 = alpha_rs[0][0];
                obj1 = model.Sum(obj1, model.Prod(-1, alpha_rs[0][0]));
                INumExpr obj2 = alpha_rs[0][0];
                obj2 = model.Sum(obj2, model.Prod(-1, alpha_rs[0][0]));
                INumExpr obj5 = alpha_rs[0][0];
                obj5 = model.Sum(obj5, model.Prod(-1, alpha_rs[0][0]));
                for (int r = 1; r < R + 1; r++)
                {
                    for (int s = 1; s < S + 1; s++)
                    {
                        obj = model.Sum(obj, model.Prod(alpha_rs[r][s], p_rs[r][s]));
                    }
                }
                for (int r = 1; r < R + 1; r++)
                {
                    for (int s = 1; s < S + 1; s++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            obj1 = model.Sum(obj1, model.Prod(theta_rsk[r][s][k], (d_ijk[0][s][k] + d_ijk[s][r][k])));
                        }
                    }
                }
                //for (int w = 0; w < W; w++)
                {
                    for (int r = 1; r < R + 1; r++)
                    {
                        for (int o = 1; o < O + 1; o++)
                        {
                            for (int k = 0; k < K; k++)
                            {
                                obj2 = model.Sum(obj2, model.Prod(d_wijk[w][r][o][k], epsilon_wrok[w][r][o][k]));
                            }
                        }
                    }
                    for (int o = 1; o < O + 1; o++)
                    {
                        for (int k = 0; k < K; k++)
                        {
                            obj5 = model.Sum(obj5, model.Prod(f, gama_wok[w][o][k]));
                        }
                    }
                }
                obj = model.Sum(obj, model.Prod(-1, obj1));
                obj2 = model.Sum(obj5, model.Prod(-1, obj2));
                obj = model.Sum(obj, obj2);
                model.AddMaximize(obj);
                #endregion

                model.SetOut(null);
                #region 输出
                if (model.Solve())
                {
                    Console.WriteLine("method3everytime = {0}   未来订单成本= {1}", Math.Round(model.GetValue(obj), 2), Math.Round(model.GetValue(obj2), 2));
                    sum[w] = Math.Round(model.GetValue(obj), 2);
                }
                else
                {
                    Console.WriteLine("无解");
                }
                #endregion

            }
            for (int w = 0; w < W; w++)
            {
                method3 = method3 + pai_w[w] * sum[w];
            }

        }
        static void Main(string[] args)
        {
            Program a = new Program();
            Random random = new Random();

            #region 参数

            string[] strArray1 = File.ReadAllLines(@"D:\wjw\WJW\实验\试验distance.txt");
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

            #region p_rs[r][s] 商家s出售已知订单r对应商品的潜在损失
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
                            d_wijk[w][i][j][k] = /*Math.Round(*/random.Next(3)/* + random.NextDouble(), 2)*/+3;
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
            double[][][] ar_ijk = new double[R + 1][][];
            double[][][] sum_ijk = new double[R + 1][][];
            for (int i = 0; i < R + 1; i++)
            {
                ar_ijk[i] = new double[O + 1][];
                sum_ijk[i] = new double[O + 1][];
                for (int j = 1; j < O + 1; j++)
                {
                    ar_ijk[i][j] = new double[K];
                    sum_ijk[i][j] = new double[K];
                    for (int k = 0; k < K; k++)
                    {
                        for (int w = 0; w < W; w++)
                        {
                            sum_ijk[i][j][k] = sum_ijk[i][j][k] + d_wijk[w][i][j][k];
                        }
                        ar_ijk[i][j][k] = Math.Round(sum_ijk[i][j][k] / W, 2);
                        //Console.WriteLine("****************sum_ijk[i{0}][j{1}][k{2}]={3}   平均{4}", i, j, k, sum_ijk[i][j][k], ar_ijk[i][j][k]);
                    }
                }
            }

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


            a.method1(p_rs, d_ijk, t_ijk, z_rs, u_r, d_wijk, pai_w);

            DateTime begin = DateTime.Now;
            double method2 = 0;
            a.method2average(ref method2, p_rs, d_ijk, t_ijk, z_rs, u_r, ar_ijk);
            a.method2second(ref method2, d_wijk, pai_w);
            Console.WriteLine("************************ method2={0}", method2);
            DateTime end = DateTime.Now;
            Console.WriteLine("######################Method2的总时间={0}", end - begin);

            double method3 = 0;
            a.method3(ref method3, p_rs, d_ijk, t_ijk, z_rs, u_r, d_wijk, pai_w);
            Console.WriteLine("************************ method3={0}", method3);

        }
    }
}
