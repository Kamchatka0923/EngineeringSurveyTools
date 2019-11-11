using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 线元法平曲线线路
    /// </summary>
    class Horizontal_LE
    {
        //连续里程
        public List<HCurve_LE> curves;

        private struct Kinfo
        {
            internal double K0;//起始里程
            internal double K1;//结束里程
            internal int index0;//起始序号
            internal int index1;//结束序号
        }
        private List<Kinfo> Krange;

        public Horizontal_LE()
        {
            curves = new List<HCurve_LE>();
            Krange = new List<Kinfo>();
        }

        /// <summary>
        /// 添加平曲线元
        /// </summary>
        /// <param name="curve">平曲线元</param>
        public void Add(HCurve_LE curve)
        {
            if (curves.Count == 0)
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = curve.K0;
                k1.index0 = 0;
                k1.K1 = curve.K1;
                k1.index1 = 0;
                Krange.Add(k1);
            }
            else if (curve.K0 == curves[curves.Count - 1].K1)
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = Krange[Krange.Count - 1].K0;
                k1.index0 = Krange[Krange.Count - 1].index0;
                k1.K1 = curve.K1;
                k1.index1 = curves.Count - 1;
                Krange.Add(k1);
            }
            else             //曲线不连续时
            {
                curves.Add(curve);
                Kinfo k1 = new Kinfo();
                k1.K0 = curve.K0;
                k1.index0 = curves.Count - 1;
                k1.K1 = curve.K1;
                k1.index1 = curves.Count - 1;
                Krange.Add(k1);
            }
        }
        /// <summary>
        /// 根据里程搜索所在曲线
        /// </summary>
        /// <param name="Kp">里程</param>
        /// <returns>曲线顺序号</returns>
        private int SearchCurveByMilleage(double Kp)
        {
            int i = 0;
            for (i = 0; i < Krange.Count; i++)
            {
                if (Kp <= Krange[i].K1 && Kp >= Krange[i].K0)//里程值在线路内
                    break;
            }
            if (i == Krange.Count)      //里程值在线路外
                return -1;
            for (i = 0; i < curves.Count; i++)
            {
                if (Kp >= curves[i].K0)
                {
                    double li = Kp - curves[i].K0;
                    if (li < curves[i].L)
                        break;
                }
            }
            return i;
        }

        /// <summary>
        /// 根据里程计算未知点坐标及方位角
        /// </summary>
        /// <param name="K">未知点里程</param>
        /// <param name="x">未知点x坐标</param>
        /// <param name="y">未知点y坐标</param>
        /// <param name="betta_p">未知点方位角</param>
        public void Milleage2XY(double K, out double x, out double y, out double betta_p)
        {
            
            int index = SearchCurveByMilleage(K);
            if (index < 0)
                throw new Exception("Horrizontal_LE:Milleage2XY:里程值不在平曲线范围内");
            switch (curves[index].type)
            {
                case 0:
                    {
                        double dx, dy;
                        double l = K - curves[index].K0;
                        dx = l * Math.Cos(curves[index].azimuth);
                        dy = l * Math.Sin(curves[index].azimuth);
                        x = curves[index].x0 + dx;
                        y = curves[index].y0 + dy;
                        betta_p = curves[index].azimuth;
                        break;
                    }
                default:
                    {
                        double xA = curves[index].x0;
                        double yA = curves[index].y0;
                        double rhoA = curves[index].rhoA;
                        double rhoB = curves[index].rhoB;
                        double L = curves[index].L;
                        double li = K - curves[index].K0;
                        double alfaA = curves[index].azimuth;
                        CoordCompute(xA, yA, rhoA, L, rhoB, li, alfaA, out x, out y);
                        betta_p = BearingCompute(rhoA, L, rhoB, li, alfaA);
                        break;
                    }
            }
        }
        /// <summary>
        /// 根据里程和偏距计算未知点平面坐标
        /// </summary>
        /// <param name="K">里程</param>
        /// <param name="offset">偏距</param>
        /// <param name="x">未知点x坐标</param>
        /// <param name="y">未知点y坐标</param>
        public void MilleageOffset2XY(double K, double offset, out double x, out double y)
        {
            Milleage2XY(K, out double xA, out double yA, out double alphaA);
            double betta = Math.PI / 2;
            if (offset != 0)
            {
                x = xA + offset * Math.Cos(alphaA + betta);
                y = yA + offset * Math.Sin(alphaA + betta);
            }
            else
            {
                x = xA;
                y = yA;
            }
        }
        /// <summary>
        /// 根据坐标搜索起点距离最近的曲线
        /// </summary>
        /// <param name="x">已知点x坐标</param>
        /// <param name="y">已知点y坐标</param>
        /// <returns></returns>
        private int SearchCurveByCoordinate(double x, double y)
        {
            int index1 = -1;
            double minDist1 = 9999999;
            double d = 0;
            double x1 = 0, y1 = 0;
            for (int i = 0; i < curves.Count; i++)
            {
                x1 = curves[i].x0;
                y1 = curves[i].y0;
                d = Dist(x1, y1, x, y);
                if (d < minDist1)
                {
                    minDist1 = d;
                    index1 = i;
                }
            }
            return index1;
        }
        /// <summary>
        /// 根据坐标计算里程和偏距
        /// </summary>
        /// <param name="x">已知点x坐标</param>
        /// <param name="y">已知点y坐标</param>
        /// <param name="K">里程</param>
        /// <param name="offset">偏距</param>
        public void XY2MilleageOffset(double x, double y,
            out double K, out double offset)
        {
            int index = SearchCurveByCoordinate(x, y);
            double xA = curves[index].x0;
            double yA = curves[index].y0;
            double alpha_A = curves[index].azimuth;
            double rhoA = curves[index].rhoA;
            double rhoB = curves[index].rhoB;
            double L = curves[index].L;
            NumberOfCompute(xA, yA, alpha_A, rhoA, rhoB, L, x, y,
                out double lp, out offset);
            if (SearchCurveByMilleage(lp + curves[index].K0 + 0.000001) < 0)
                throw new Exception("Horizontal_LE:XY2MilleageOffset:该坐标位于平曲线外！");
            if ((lp < 0) && (index > 0))
            {
                index = index - 1;
                xA = curves[index].x0;
                yA = curves[index].y0;
                alpha_A = curves[index].azimuth;
                rhoA = curves[index].rhoA;
                rhoB = curves[index].rhoB;
                L = curves[index].L;
                NumberOfCompute(xA, yA, alpha_A, rhoA, rhoB, L, x, y,
                out lp, out offset);
            }
            K = lp + curves[index].K0;
        }


        /// <summary>
        /// 计算回旋曲线上任意一点的切线方位角(弧度制)
        /// </summary>
        /// <param name="rhoA">回旋曲线起点曲率</param>
        /// <param name="L">起点终点里程差</param>
        /// <param name="rhoB">回旋曲线终点曲率</param>
        /// <param name="li">待求点到回旋曲线起点的长度</param>
        /// <param name="alfaA">起点处切线方位角(弧度制)</param>
        /// <returns></returns>
        private double BearingCompute(double rhoA, double L, double rhoB,
            double li, double alfaA)
        {
            double rhoi;        //回旋曲线上任一点的曲率
            rhoi = rhoA + (rhoB - rhoA) * li / L;
            double alfa_i;      //回旋曲线上一点处的切线方位角
            double betta_i;     //曲线上一点处的回旋角

            betta_i = (rhoA + rhoi) * li / 2;
            alfa_i = alfaA + rhoA * li + (rhoB - rhoA) * li * li / (2 * L);
            return alfa_i;
        }
        /// <summary>
        /// 计算曲线上任意一点的坐标（迭代）
        /// </summary>
        /// <param name="xA">回旋曲线起点x坐标</param>
        /// <param name="yA">回旋曲线起点y坐标</param>
        /// <param name="rhoA">回旋曲线起点曲率</param>
        /// <param name="L">起点终点里程差</param>
        /// <param name="rhoB">回旋曲线终点曲率</param>
        /// <param name="li">待求点到回旋曲线起点的长度</param>
        /// <param name="alfaA">起点处切线方位角(弧度制)</param>
        /// <param name="xi">回旋曲线上任一点的x坐标</param>
        /// <param name="yi">回旋曲线上任一点的y坐标</param>
        private void CoordCompute(double xA, double yA,
            double rhoA, double L, double rhoB, double li, double alfaA,
            out double xi, out double yi)
        {
            xi = 0; yi = 0;
            double alfa_i, alfa_k;
            alfa_i = BearingCompute(rhoA, L, rhoB, li, alfaA);
            int m, n, k;
            double xi_1, yi_1;  //定义坐标的暂存变量
            m = 1;
            double temp;
            do
            {
                xi_1 = xi;
                yi_1 = yi;
                n = m * 2;
                double temp1, temp2, temp3, temp4;  //计算坐标的暂存变量
                temp1 = 0;
                temp2 = 0;
                temp3 = 0;
                temp4 = 0;
                double h, lk;
                h = li / n;
                for (k = 1; k <= m; k++)
                {
                    lk = h * (2 * k - 1);
                    alfa_k = BearingCompute(rhoA, L, rhoB, lk, alfaA);
                    temp1 = temp1 + Math.Cos(alfa_k);
                    temp3 = temp3 + Math.Sin(alfa_k);
                }
                for (int j = 1; j <= m - 1; j++)
                {
                    lk = h * (2 * j);
                    alfa_k = BearingCompute(rhoA, L, rhoB, lk, alfaA);
                    temp2 = temp2 + Math.Cos(alfa_k);
                    temp4 = temp4 + Math.Sin(alfa_k);
                }
                xi = xA + h / 3 * (Math.Cos(alfaA) + Math.Cos(alfa_i) +
                    4 * temp1 + 2 * temp2);
                yi = yA + h / 3 * (Math.Sin(alfaA) + Math.Sin(alfa_i) +
                    4 * temp3 + 2 * temp4);
                temp = Math.Sqrt((xi - xi_1) * (xi - xi_1) + (yi - yi_1) * (yi - yi_1));
                m++;
            } while (temp > 0.00001);
        }
        /// <summary>
        /// 由曲线外一点反求桩号
        /// </summary>
        /// <param name="xA">曲线元起点A的x坐标</param>
        /// <param name="yA">曲线元起点A的y坐标</param>
        /// <param name="alfa_A">起点切线方位角（弧度制）</param>
        /// <param name="rho_A">起点曲率</param>
        /// <param name="rho_B">终点曲率</param>
        /// <param name="L">曲线元长度</param>
        /// <param name="xp1">路线外P点的x坐标</param>
        /// <param name="yp1">路线外P点的y坐标</param>
        /// <param name="lp">P点距起点的里程差</param>
        /// <param name="Dp">P点到中桩的距离，Dp小于0时，P点位于线路左侧，Dp大于0时，P点位于线路右侧</param>
        private void NumberOfCompute(double xA, double yA, double alfa_A,
            double rho_A, double rho_B, double L, double xp1, double yp1,
            out double lp, out double Dp)
        {
            double xi, yi, alfa_i;
            xi = xA;
            yi = yA;
            alfa_i = alfa_A;
            lp = 0;
            double dA;      //迭代中产生的距离
            double rhoB;    //迭代中产生的方位角
            do
            {
                dA = (yp1 - yi) * Math.Cos(alfa_i - Math.PI / 2) -
                    (xp1 - xi) * Math.Sin(alfa_i - Math.PI / 2);
                lp = dA + lp;
                alfa_i = BearingCompute(rho_A, L, rho_B, lp, alfa_A);
                CoordCompute(xA, yA, rho_A, L, rho_B, lp, alfa_A, out xi, out yi);
            } while (dA > 0.001);
            Dp = (xi - xp1) / Math.Sin(alfa_i);
        }
        /// <summary>
        /// 求平面上两点距离
        /// </summary>
        private double Dist(double X1, double Y1, double X2, double Y2)
        {
            double d;
            d = Math.Sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
            return d;
        }
    }
}
