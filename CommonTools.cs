using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SurveyTools
{
    /// <summary>
    /// 角度转换
    /// 交会计算
    /// 测量工具
    /// 线路平纵曲线计算
    /// 坐标转换
    /// </summary>
    public class CommonTools
    {
        #region 角度转换：度分秒、十进制、弧度制之间的转换
        /// <summary>
        /// 度分秒转十进制度
        /// </summary>
        /// <param name="dDms">度分秒格式角度</param>
        /// <returns>十进制格式角度</returns>
        public static double Dms2Deg(double dDms)
        {
            int iDegree, iMin;
            double dSec, dDeg;
            iDegree = (int)dDms;
            iMin = (int)((dDms - iDegree) * 100);
            dSec = ((dDms - iDegree) * 100 - iMin) * 100;
            dDeg = iDegree + ((double)iMin / 60) + dSec / 3600;
            return dDeg;
        }
        /// <summary>
        /// 十进制度转度分秒
        /// </summary>
        /// <param name="dDeg">十进制格式角度</param>
        /// <returns>度分秒格式角度</returns>
        public static double Deg2Dms(double dDeg)
        {
            int iDegree, iMin;
            double dSec, dTmp, dDms;
            iDegree = (int)dDeg;
            dTmp = (dDeg - iDegree) * 60;
            iMin = (int)dTmp;
            dSec = (dTmp - iMin) * 60;
            dDms = iDegree + (double)iMin / 100 + dSec / 10000;
            return dDms;
        }
        /// <summary>
        /// 十进制度转弧度制
        /// </summary>
        /// <param name="dDeg">十进制格式角度</param>
        /// <returns>弧度制</returns>
        public static double Deg2Rad(double dDeg)
        {
            double dRad;
            dRad = (dDeg / 180) * Math.PI;
            return dRad;
        }
        /// <summary>
        /// 弧度制转十进制度
        /// </summary>
        /// <param name="dRad">弧度制</param>
        /// <returns>十进制格式角度</returns>
        public static double Rad2Deg(double dRad)
        {
            double dDeg;
            dDeg = 180 * (dRad / Math.PI);
            return dDeg;
        }
        /// <summary>
        /// 度分秒转弧度制
        /// </summary>
        /// <param name="dDms">度分秒格式角度</param>
        /// <returns>弧度制</returns>
        public static double Dms2Rad(double dDms)
        {
            double dRad, dSec;
            int iDegree, iMin;
            iDegree = (int)dDms;
            iMin = (int)((dDms - iDegree) * 100);
            dSec = ((dDms - iDegree) * 100 - iMin) * 100;
            dRad = Math.PI * (iDegree + ((double)iMin / 60) + dSec / 3600) / 180;
            return dRad;
        }
        /// <summary>
        /// 弧度制转度分秒
        /// </summary>
        /// <param name="dRad">弧度制</param>
        /// <returns>度分秒格式角度</returns>
        public static double Rad2Dms(double dRad)
        {
            int iDegree, iMin;
            double dSec, dDms;
            iDegree = (int)((dRad / Math.PI) * 180);
            iMin = (int)(((dRad / Math.PI) * 180 - iDegree) * 60);
            dSec = (((dRad / Math.PI) * 180 - iDegree) * 60 - iMin) * 60;
            dDms = iDegree + (double)iMin / 100 + dSec / 10000;
            return dDms;
        }
        #endregion

        #region 交会计算：前方、后方、边长交会
        /// <summary>
        /// 前方交会计算（A、B、P逆时针）
        /// </summary>
        /// <param name="Xa">已知点a坐标</param>
        /// <param name="Ya">已知点a坐标</param>
        /// <param name="Xb">已知点b坐标</param>
        /// <param name="Yb">已知点b坐标</param>
        /// <param name="Alfa">a点交会角（弧度）</param>
        /// <param name="Beta">b点交会角（弧度）</param>
        /// <param name="Xp">待定点p坐标</param>
        /// <param name="Yp">待定点p坐标</param>
        public static void ForeIntersecPos(double Xa, double Ya, double Xb, double Yb,
            double Alfa, double Beta, out double Xp, out double Yp)
        {
            double ctgA, ctgB;
            //反正切
            ctgA = 1 / Math.Tan(Alfa);
            ctgB = 1 / Math.Tan(Beta);
            //计算前方交会定位值
            Xp = ((Xa * ctgB + Xb * ctgA) + (Yb - Ya)) / (ctgA + ctgB);
            Yp = ((Ya * ctgB + Yb * ctgA) + (Xa - Xb)) / (ctgA + ctgB);
        }
        /// <summary>
        /// 测角后方交会
        /// </summary>
        /// <param name="Xa"></param>
        /// <param name="Ya"></param>
        /// <param name="Xb"></param>
        /// <param name="Yb"></param>
        /// <param name="Xc"></param>
        /// <param name="Yc"></param>
        /// <param name="alpha">BC边对应测角</param>
        /// <param name="betta">AC边对应测角</param>
        /// <param name="gamma">AB边对应测角</param>
        /// <param name="Xp"></param>
        /// <param name="Yp"></param>
        public static void ResIntersecPos(double Xa, double Ya, double Xb, double Yb,
            double Xc, double Yc, double alpha, double betta, double gamma,
            out double Xp, out double Yp)
        {
            double A, B, C, a, b, c;
            a = Dist(Xc, Yc, Xb, Yb);
            b = Dist(Xc, Yc, Xa, Ya);
            c = Dist(Xa, Ya, Xb, Yb);
            A = Math.Acos((b * b + c * c - a * a) / (2 * b * c));
            B = Math.Acos((a * a + c * c - b * b) / (2 * a * c));
            C = Math.Acos((b * b + a * a - c * c) / (2 * b * a));
            double c_temp, b_temp, a_temp;//判断是否在危险圆附近的判断变量
            c_temp = Rad2Deg(alpha + betta + C);
            b_temp = Rad2Deg(alpha + gamma + B);
            a_temp = Rad2Deg(betta + gamma + A);
            if ((c_temp > 170 && c_temp < 190) || (c_temp > 170 && c_temp < 190) ||
                (c_temp > 170 && c_temp < 190))
                throw new Exception("Tools:ResIntersecPos:在危险圆附近");
            double pa = (Math.Tan(alpha) * Math.Tan(A)) / (Math.Tan(alpha) - Math.Tan(A));
            double pb = (Math.Tan(betta) * Math.Tan(B)) / (Math.Tan(betta) - Math.Tan(B));
            double pc = (Math.Tan(gamma) * Math.Tan(C)) / (Math.Tan(gamma) - Math.Tan(C));
            Xp = (pa * Xa + pb * Xb + pc * Xc) / (pa + pb + pc);
            Yp = (pa * Ya + pb * Yb + pc * Yc) / (pa + pb + pc);
        }
        /// <summary>
        /// 边长交会
        /// </summary>
        /// <param name="Xa"></param>
        /// <param name="Ya"></param>
        /// <param name="Xb"></param>
        /// <param name="Yb"></param>
        /// <param name="Dap">AP距离</param>
        /// <param name="Dbp">BP距离</param>
        /// <param name="Xp"></param>
        /// <param name="Yp"></param>
        public static void LinearIntersecPos(double Xa, double Ya, double Xb, double Yb,
            double Dap, double Dbp, out double Xp, out double Yp)
        {
            double Dab = Dist(Xa, Ya, Xb, Yb);
            double alpha_ab = Azimuth(Xa, Ya, Xb, Yb);
            double A = Math.Acos((Dab * Dab + Dap * Dap - Dbp * Dbp) / (2 * Dab * Dap));
            double alpha_ap = alpha_ab - A;
            Xp = Xa + Dap * Math.Cos(alpha_ap);
            Yp = Ya + Dap * Math.Sin(alpha_ap);
        }
        #endregion

        #region 测量工具
        /// <summary>
        /// 判断一个平面点在测量坐标系中的象限位置
        /// </summary>
        /// <param name="GeoX">平面点横坐标</param>
        /// <param name="GeoY">平面点纵坐标</param>
        /// <returns>
        /// 0：原点
        /// 1：第一象限
        /// 2：第一象限
        /// 3：第一象限
        /// 4：第一象限
        /// 104：X正半轴
        /// 102：X负半轴
        /// 203：Y正半轴
        /// 304：Y负半轴
        /// </returns>
        public static int PointJudgeQuadrant(double GeoX, double GeoY)
        {
            int nQuadrant = -1;
            if (GeoX > 0)
            {
                if (GeoY > 0)
                    nQuadrant = 1;      //第一象限
                else if (GeoY < 0)
                    nQuadrant = 4;      //第四象限
                else
                    nQuadrant = 104;    //x正半轴
            }
            else if (GeoX < 0)
            {
                if (GeoY > 0)
                    nQuadrant = 2;      //第二象限
                else if (GeoY < 0)
                    nQuadrant = 3;      //第三象限
                else
                    nQuadrant = 203;    //x负半轴
            }
            else if (GeoY > 0)
                nQuadrant = 103;        //y正半轴
            else if (GeoY > 0)
                nQuadrant = 304;        //y负半轴
            else if (GeoY > 0)
                nQuadrant = 0;          //原点
            return nQuadrant;
        }
        /// <summary>
        /// 同一参考椭球下的三维地心坐标（笛卡尔坐标系）转换为大地坐标
        /// 东经0到180，Y大于0;西经0到180，Y小于0 
        /// </summary>
        /// <param name="X">大地坐标X</param>
        /// <param name="Y">大地坐标Y</param>
        /// <param name="Z">大地坐标Z</param>
        /// <param name="a">椭球长半径</param>
        /// <param name="e">椭球偏心率</param>
        /// <param name="dB">大地纬度（弧度制）</param>
        /// <param name="dL">大地经度（弧度制）</param>
        /// <param name="dH">大地高</param>
        public static void Descartes2Geodetic(double X, double Y, double Z,
            double a, double e, out double dB, out double dL, out double dH)
        {
            dL = Math.Atan(Y / X);
            double dH0 = 0;
            double dS = Math.Sqrt(X * X + Y * Y);
            double dB0 = Math.Atan(Z / (dS * (1 - e * e)));
            double dh;
            do
            {
                double W = Math.Sqrt(1 - e * e * Math.Sin(dB0) * Math.Sin(dB0));
                double N = a / W;
                dH = dS / Math.Cos(dB0) - N;    //计算大地高新值
                dB = Math.Atan(Z / (dS * (1 - e * e * N / (N + dH))));//计算纬度新值
                dh = dH - dH0;  //计算大地高迭代新旧值之差

                //新旧值替换
                dH0 = dH;
                dB0 = dB;
            } while (Math.Abs(dh) >= 0.0001);
        }
        /// <summary>
        /// 方位角计算
        /// </summary>
        /// <param name="X1">已知点1的X坐标</param>
        /// <param name="Y1">已知点1的Y坐标</param>
        /// <param name="X2">已知点2的X坐标</param>
        /// <param name="Y2">已知点2的Y坐标</param>
        /// <returns>方位角（弧度制）</returns>
        public static double Azimuth(double X1, double Y1, double X2, double Y2)
        {
            int sgn = 0;
            double dx, dy;
            double EPSILON = 1.0E-10;
            dx = X2 - X1;
            dy = Y2 - Y1 + EPSILON;
            if (dy >= 0)
                sgn = 1;
            else
                sgn = -1;
            return Math.PI - sgn * Math.PI / 2 - Math.Atan(dx / dy);
        }
        /// <summary>
        /// 求平面上两点距离
        /// </summary>
        public static double Dist(double X1, double Y1, double X2, double Y2)
        {
            double d;
            d = Math.Sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
            return d;
        }
        /// <summary>
        /// 求空间上两点距离
        /// </summary>
        public static double Dist(double X1, double Y1, double Z1,
            double X2, double Y2, double Z2)
        {
            double d;
            d = Math.Sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) +
                (Z2 - Z1) * (Z2 - Z1));
            return d;
        }
        /// <summary>
        /// 求三角形的三个内角
        /// </summary>
        /// <param name="Xa">顶点a的X坐标</param>
        /// <param name="Ya">顶点a的Y坐标</param>
        /// <param name="Xb">顶点b的X坐标</param>
        /// <param name="Yb">顶点b的Y坐标</param>
        /// <param name="Xc">顶点c的X坐标</param>
        /// <param name="Yc">顶点c的Y坐标</param>
        /// <param name="a">a角度（弧度）</param>
        /// <param name="b">b角度（弧度）</param>
        /// <param name="c">c角度（弧度）</param>
        public static void GetInnerAngle(double Xa, double Ya, double Xb, double Yb,
            double Xc, double Yc, out double a, out double b, out double c)
        {
            double Sa, Sb, Sc, cosa, cosb, cosc;
            Sa = Dist(Xc, Yc, Xb, Yb);
            Sb = Dist(Xc, Yc, Xa, Ya);
            Sc = Dist(Xa, Ya, Xb, Yb);
            cosa = (Sb * Sb + Sc * Sc - Sa * Sa) / (2 * (Sb * Sc));
            cosb = (Sa * Sa + Sc * Sc - Sb * Sb) / (2 * (Sa * Sc));
            cosc = (Sa * Sa + Sb * Sb - Sc * Sc) / (2 * (Sa * Sb));
            a = Math.Acos(cosa);
            b = Math.Acos(cosb);
            c = Math.Acos(cosc);
        }
        /// <summary>
        /// 计算三角形外接圆的圆心坐标
        /// </summary>
        /// <param name="Xa"></param>
        /// <param name="Ya"></param>
        /// <param name="Xb"></param>
        /// <param name="Yb"></param>
        /// <param name="Xc"></param>
        /// <param name="Yc"></param>
        /// <param name="X0">圆心X坐标</param>
        /// <param name="Y0">圆心Y坐标</param>
        public static void CycleCenter(double Xa, double Ya, double Xb, double Yb,
            double Xc, double Yc, out double X0, out double Y0)
        {
            double a, b, c, g;
            a = Xa * Xa + Ya * Ya;
            b = Xb * Xb + Yb * Yb;
            c = Xc * Xc + Yc * Yc;
            g = (Yc - Yb) * Xa + (Ya - Yc) * Xb + (Yb - Ya) * Xc;
            X0 = ((b - c) * Ya + (c - a) * Yb + (a - b) * Yc) / (2 * g);
            Y0 = ((b - c) * Xa + (c - a) * Xb + (a - b) * Xc) / (2 * g);
        }
        /// <summary>
        /// 边长距离改化
        /// </summary>
        /// <param name="y1">椭球体上P1点横坐标</param>
        /// <param name="y2">椭球体上P2点横坐标</param>
        /// <param name="S">大地线长度</param>
        /// <param name="Rm">大地线始末端点的平均纬度计算（查取）的椭球平均曲率半径</param>
        /// <returns></returns>
        public static double Cdistance_correction(double y1, double y2, double S, double Rm)
        {
            double ym;
            double dy;
            ym = (y1 + y2) / 2;
            dy = y1 - y2;
            double D = (1 + (ym * ym) / (2 * Rm * Rm) + (dy * dy) / (24 + Rm + Rm) + Math.Pow(ym, 4) /
                (24 * Math.Pow(Rm, 4))) * S;
            return D;
        }
        /// <summary>
        /// 三角高程测量计算
        /// </summary>
        /// <param name="D">水平距离</param>
        /// <param name="Ha">站点高程</param>
        /// <param name="i">仪器高</param>
        /// <param name="v">目标高</param>
        /// <param name="alpha">观测垂直角（度分秒）</param>
        /// <param name="Hb">目标点高程</param>
        public static void Ctrigonometric_leveling(double D, double Ha, double i,
            double v, double alpha, out double Hb)
        {
            if (D < 0 || i < 0 || v < 0)
                throw new Exception("Tools:Ctrigonometric_leveling:数据格式错误");
            double k = 0.14;    //大气折光系数
            double R = 6371000; //地球平均曲率半径
            double h_ab;        //A、B两点的高差
            alpha = Dms2Rad(alpha);
            double f;           //两差改正
            f = (1 - k) * D * D / (2.0 * R);
            h_ab = D * Math.Tan(alpha) + i - v + f;
            Hb = Ha + h_ab;
        }
        #endregion

        #region 线路平纵计算工具
        /// <summary>
        /// 完整回旋曲线上任意一点的坐标计算
        /// </summary>
        /// <param name="R">回旋曲线终点处曲率半径</param>
        /// <param name="ls">回旋曲线总长度</param>
        /// <param name="l">回旋曲线上点距离起点的长度</param>
        /// <param name="xp">该点x坐标</param>
        /// <param name="yp">该点y坐标</param>
        public static void CoordOfPointInCurve(double R, double ls, double l,
            out double xp, out double yp)
        {
            xp = l - Math.Pow(l, 5) / (40 * R * R * ls * ls) +
                Math.Pow(l, 9) / (3456 * Math.Pow(R, 4) * Math.Pow(ls, 4));
            yp = Math.Pow(l, 3) / (6 * R * ls) -
                Math.Pow(l, 7) / (336 * Math.Pow(R, 3) * Math.Pow(ls, 3) +
                Math.Pow(l, 11) / (42240 * Math.Pow(R, 5)) * Math.Pow(ls, 5));
        }
        /// <summary>
        /// 完整回旋曲线终点坐标计算（由终点处曲率半径计算）
        /// </summary>
        /// <param name="R">回旋曲线终点处曲率半径</param>
        /// <param name="ls">回旋曲线总长度</param>
        /// <param name="xp">终点x坐标</param>
        /// <param name="yp">终点y坐标</param>
        public static void EndPointCoord_RAndls(double R, double ls,
            out double xp, out double yp)
        {
            xp = ls - Math.Pow(ls, 3) / (40 * R * R) +
                Math.Pow(ls, 5) / (3456 * Math.Pow(R, 4));
            yp = Math.Pow(ls, 2) / (6 * R) - Math.Pow(ls, 4) / (336 * Math.Pow(R, 3)) +
                Math.Pow(ls, 6) / (42240 * Math.Pow(R, 5));
        }
        /// <summary>
        /// 完整回旋曲线终点坐标计算（由回旋曲线参数计算）
        /// </summary>
        /// <param name="A">回旋曲线参数</param>
        /// <param name="ls">回旋曲线总长度</param>
        /// <param name="xp">终点x坐标</param>
        /// <param name="yp">终点y坐标</param>
        public static void EndPointCoord_AAndls(double A, double ls,
            out double xp, out double yp)
        {
            xp = ls - Math.Pow(ls, 5) / (40 * Math.Pow(A, 4)) +
                Math.Pow(ls, 9) / (3456 * Math.Pow(A, 8));
            yp = Math.Pow(ls, 3) / (6 * A * A) - Math.Pow(ls, 7) / (336 * Math.Pow(A, 6)) +
                Math.Pow(ls, 11) / (42240 * Math.Pow(A, 10));
        }
        /// <summary>
        /// 线路基本型曲线要素计算
        /// </summary>
        /// <param name="ls1">前回旋曲线长度</param>
        /// <param name="ls2">后回旋曲线长度</param>
        /// <param name="R">圆曲线半径</param>
        /// <param name="alpha">转角（度分秒）</param>
        /// <param name="p1">前回旋曲线内移值</param>
        /// <param name="p2">后回旋曲线内移值</param>
        /// <param name="q1">前回旋曲线切线增长值</param>
        /// <param name="q2">后回旋曲线切线增长值</param>
        /// <param name="T1">前回旋曲线切线长度</param>
        /// <param name="T2">后回旋曲线切线长度</param>
        /// <param name="Ly">圆曲线长</param>
        /// <param name="L">平曲线长</param>
        /// <param name="E">外距</param>
        public static void ParaOfRoad(double ls1, double ls2, double R, double alpha,
            out double p1, out double p2, out double q1, out double q2,
            out double T1, out double T2, out double Ly, out double L, out double E)
        {
            alpha = Math.Abs(alpha);
            alpha = Dms2Rad(alpha);
            if ((R * alpha) == 0)
                throw new Exception("Tools:ParaOfRoad:请输入正确的参数！");
            else
            {
                p1 = ls1 * ls1 / (24 * R) - Math.Pow(ls1, 4) / (2688 * Math.Pow(R, 3)) +
                    Math.Pow(ls1, 6) / (506880 * Math.Pow(R, 5));
                p2 = ls2 * ls2 / (24 * R) - Math.Pow(ls2, 4) / (2688 * Math.Pow(R, 3)) +
                    Math.Pow(ls2, 6) / (506880 * Math.Pow(R, 5));
                q1 = ls1 / 2 - Math.Pow(ls1, 3) / (240 * R * R) +
                    Math.Pow(ls1, 5) / (34560 * Math.Pow(R, 4));
                q2 = ls2 / 2 - Math.Pow(ls2, 3) / (240 * R * R) +
                    Math.Pow(ls2, 5) / (34560 * Math.Pow(R, 4));
                T1 = q1 + (R + p2 - (R + p1) * Math.Cos(alpha)) / Math.Sin(alpha);
                T2 = q2 + (R + p1 - (R + p2) * Math.Cos(alpha)) / Math.Sin(alpha);
                Ly = R * alpha - ls1 / 2 - ls2 / 2;
                L = Ly + ls1 + ls2;
                E = (R + (p1 + p2) / 2) * (1 / Math.Cos(alpha / 2)) - R;
            }
        }
        /// <summary>
        /// 曲线主点桩号计算
        /// </summary>
        /// <param name="Kjd">交点里程</param>
        /// <param name="T1">前回旋曲线切线长度</param>
        /// <param name="ls1">前回旋曲线长度</param>
        /// <param name="Ly">圆曲线长</param>
        /// <param name="ls2">后回旋曲线长度</param>
        /// <param name="Kzh">前回旋曲线起点里程</param>
        /// <param name="Khy">前回旋曲线终点里程</param>
        /// <param name="Kqz">圆曲线中点里程</param>
        /// <param name="Kyh">后回旋曲线起点里程</param>
        /// <param name="Khz">后回旋曲线终点里程</param>
        public static void NumberOfMainPoint(double Kjd, double T1, double ls1,
            double Ly, double ls2,
            out double Kzh, out double Khy, out double Kqz, out double Kyh, out double Khz)
        {
            Kzh = Kjd - T1;
            Khy = Kzh + ls1;
            Kqz = Khy + Ly / 2;
            Kyh = Khy + Ly;
            Khz = Kyh + ls2;
        }
        /// <summary>
        /// 线路基本型中线坐标计算
        /// </summary>
        /// <param name="xi">平面交点x坐标</param>
        /// <param name="yi">平面交点y坐标</param>
        /// <param name="R">圆曲线半径</param>
        /// <param name="ls1">前回旋曲线长度</param>
        /// <param name="ls2">后回旋曲线长度</param>
        /// <param name="T1">前切线长度</param>
        /// <param name="T2">后切线长度</param>
        /// <param name="Ai_1">前直线方位角(弧度制)</param>
        /// <param name="Ai_2">后直线方位角(弧度制)</param>
        /// <param name="Kp">待求点中桩里程</param>
        /// <param name="Kzh">前回旋曲线起点里程（ZH点里程）</param>
        /// <param name="Khz">后回旋曲线终点里程（HZ点里程）</param>
        /// <param name="Khy">前回旋曲线终点里程（HY点里程）</param>
        /// <param name="Kyh">后回旋曲线起点里程（YH点里程）</param>
        /// <param name="alpha">转角（度分秒）</param>
        /// <param name="xp">待求点x坐标</param>
        /// <param name="yp">待求点y坐标</param>
        /// <param name="betta_p">切线方位角</param>
        public static void CoordAndBearing(double xi, double yi, double R,
            double ls1, double ls2, double T1, double T2, double Ai_1, double Ai_2,
            double Kp, double Kzh, double Khz, double Khy, double Kyh, double alpha,
            out double xp, out double yp, out double betta_p)
        {
            double l;
            //回旋曲线与直线段交点的坐标
            double Xzh, Yzh;
            double Xhz, Yhz;
            double I;
            if (alpha <= 0)
                I = -1;
            else
                I = 1;
            xp = 0;
            yp = 0;
            betta_p = 0;
            if (Kp <= Kzh)
            {
                l = Kzh - Kp;
                xp = xi - (T1 + l) * Math.Cos(Ai_1);
                yp = yi - (T1 + l) * Math.Sin(Ai_1);
                betta_p = Rad2Dms(Ai_1);
            }
            else if (Kp >= Khz)
            {
                l = Kp - Khz;
                xp = xi - (T2 + l) * Math.Cos(Ai_2);
                yp = yi - (T2 + l) * Math.Sin(Ai_2);
                betta_p = Rad2Dms(Ai_2);
            }
            else if ((Kp >= Kzh && Kp <= Khy) || (Kp >= Khy && Kp <= Kyh))
            //第一回旋曲线和圆曲线内的坐标和方位角
            {
                Xzh = xi - T1 * Math.Cos(Ai_1);
                Yzh = yi - T1 * Math.Sin(Ai_1);
                double xp_1 = 0, yp_1 = 0;      //P点局部坐标
                l = Kp - Kzh;
                if (Kp >= Kzh && Kp <= Khy)    //第一回旋曲线内的坐标和方位角
                {
                    xp_1 = l - Math.Pow(l, 5) / (40 * R * R * ls1 * ls1) +
                        Math.Pow(l, 9) / (3456 * Math.Pow(R, 4) * Math.Pow(ls1, 4));
                    yp_1 = Math.Pow(l, 3) / (6 * R * ls1) -
                        Math.Pow(l, 7) / (336 * Math.Pow(R, 3) * Math.Pow(ls1, 3)) +
                        Math.Pow(l, 11) / (42240 * Math.Pow(R, 5) * Math.Pow(ls1, 5));
                    betta_p = Ai_1 + I * l * l / (2 * ls1 * R);
                    betta_p = Rad2Dms(betta_p);
                }
                else if (Kp >= Khy && Kp <= Kyh)    //圆曲线内的坐标和方位角
                {
                    double betta;
                    betta = (2 * l - ls1) / (2 * R);
                    double p1, q1;
                    p1 = ls1 * ls1 / (24 * R) - Math.Pow(ls1, 4) / (2688 * Math.Pow(R, 3)) +
                        Math.Pow(ls1, 6) / (506880 * Math.Pow(R, 5));
                    q1 = ls1 / 2 - Math.Pow(ls1, 3) / (240 * R * R) +
                        Math.Pow(ls1, 5) / (34560 * Math.Pow(R, 4));
                    xp_1 = R * Math.Sin(betta) + q1;
                    yp_1 = R * (l - Math.Cos(betta)) + p1;
                    betta_p = Rad2Dms(Ai_1 + I * betta);
                }
                xp = Xzh + xp_1 * Math.Cos(Ai_1) - I * yp_1 * Math.Sin(Ai_1);
                yp = Yzh + xp_1 * Math.Sin(Ai_1) + I * yp_1 * Math.Cos(Ai_1);
            }
            else if (Kp >= Kyh && Kp <= Khz)
            {
                Xhz = xi + T2 * Math.Cos(Ai_2);
                Yhz = yi + T2 * Math.Sin(Ai_2);
                l = Khz - Kp;
                double xp_2, yp_2;  //P点局部坐标
                xp_2 = l - Math.Pow(l, 5) / (40 * R * R * ls2 * ls2) +
                        Math.Pow(l, 9) / (3456 * Math.Pow(R, 4) * Math.Pow(ls2, 4));
                yp_2 = Math.Pow(l, 3) / (6 * R * ls2) -
                    Math.Pow(l, 7) / (336 * Math.Pow(R, 3) * Math.Pow(ls2, 3)) +
                    Math.Pow(l, 11) / (42240 * Math.Pow(R, 5) * Math.Pow(ls2, 5));
                xp = Xhz - xp_2 * Math.Cos(Ai_2) - I * yp_2 * Math.Sin(Ai_2);
                yp = Yhz - xp_2 * Math.Sin(Ai_2) + I * yp_2 * Math.Cos(Ai_2);
                betta_p = Ai_2 - I * l * l / (2 * ls2 * R);
                betta_p = Rad2Dms(betta_p);
            }
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
        public static double BearingCompute(double rhoA, double L, double rhoB,
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
        /// 计算回旋曲线上任意一点的切线方位角(弧度制)
        /// </summary>
        /// <param name="rhoA">回旋曲线起点曲率</param>
        /// <param name="rhoB">回旋曲线终点曲率</param>
        /// <param name="KA">回旋曲线起点里程</param>
        /// <param name="KB">回旋曲线终点里程</param>
        /// <param name="Ki">待求点里程</param>
        /// <param name="alfaA">起点处切线方位角(弧度制)</param>
        /// <returns></returns>
        public static double BearingCompute(double rhoA, double rhoB,
            double KA, double KB, double Ki, double alfaA)
        {
            double L, li;
            L = KB - KA;        //线路总里程差
            li = Ki - KA;       //任一点距回旋曲线起点的长度
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
        public static void CoordCompute(double xA, double yA,
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
        /// 计算曲线上任意一点的坐标（不迭代）
        /// </summary>
        /// <param name="xA">回旋曲线起点x坐标</param>
        /// <param name="yA">回旋曲线起点y坐标</param>
        /// <param name="RA">起点曲率半径，最大为double.MaxValue</param>
        /// <param name="RB">终点曲率半径，最大为double.MaxValue</param>
        /// <param name="L">起点终点里程差</param>
        /// <param name="l">待求点到起点的里程差</param>
        /// <param name="alpha">起点处切线方位角（弧度制）</param>
        /// <param name="x">回旋曲线上任一点的x坐标</param>
        /// <param name="y">回旋曲线上任一点的y坐标</param>
        public static void CoordCompute1(double xA, double yA, double RA, double RB,
            double L, double l,double alpha, out double x, out double y)
        {
            double A = 1 / RA;
            double B = (1 / RB - 1 / RA) / (2 * L);
            y = A * l * l / 2 + B * l * l * l / 3 - Math.Pow(A, 3) * Math.Pow(l, 4) / 24 -
                A * A * B * Math.Pow(l, 5) / 10 -
                (A * B * B / 12 - Math.Pow(A, 5) / 720) * Math.Pow(l, 6) -
                (B * B * B / 42 - Math.Pow(A, 4) * B / 168) * Math.Pow(l, 7) +
                A * A * A * B * B * Math.Pow(l, 8) / 96 +
                A * A * B * B * B * Math.Pow(l, 9) / 108;
            x = l - A * A * l * l * l / 6 - A * B * Math.Pow(l, 4) / 4 -
                (B * B / 10 - Math.Pow(A, 4) / 120) * Math.Pow(l, 5) +
                A * A * A * B * Math.Pow(l, 6) / 36 +
                A * A * B * B * Math.Pow(l, 7) / 28 +
                A * B * B * B * Math.Pow(l, 8) / 48 +
                Math.Pow(B, 4) * Math.Pow(l, 9) / 216;
            CommonTools.HPOINT SourcePoint = new CommonTools.HPOINT();
            SourcePoint.X = x;
            SourcePoint.Y = y;
            CommonTools.FourPara para = new CommonTools.FourPara();
            para.DeltaX = xA;
            para.DeltaY = yA;
            para.Scale = 1;
            para.Sita = alpha;
            CommonTools.HPOINT p1 = CommonTools.HPointTrans(SourcePoint, para);
            x = p1.X;
            y = p1.Y;
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
        public static void NumberOfCompute(double xA, double yA, double alfa_A,
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
        /// 计算竖曲线要素
        /// </summary>
        /// <param name="i1">后坡坡度（弧度制）</param>
        /// <param name="i2">前坡坡度（弧度制）</param>
        /// <param name="R">竖曲线半径</param>
        /// <param name="K">变坡点桩号</param>
        /// <param name="H">变坡点高程</param>
        /// <param name="L">竖曲线长度</param>
        /// <param name="T">切线长度</param>
        /// <param name="E">竖曲线外距</param>
        /// <param name="Szy">起点桩号</param>
        /// <param name="Syz">终点桩号</param>
        public static void ParaOfVerticurve(double i1, double i2, double R, double K, double H,
            out double L, out double T, out double E, out double Szy, out double Syz)
        {
            double omi;
            if (R == 0)
                throw new Exception("Tools:ParaOfVerticurve:请输入正确的参数");
            else
            {
                omi = Math.Abs(i1 - i2);
                L = R * omi;
                T = L / 2;
                E = T * T / (2 * R);
                Szy = K - T;
                Syz = K + T;
            }
        }
        /// <summary>
        /// 计算竖曲线上一点的高程
        /// </summary>
        /// <param name="Kp">该点桩号</param>
        /// <param name="K">变坡点桩号</param>
        /// <param name="Szy">起点桩号</param>
        /// <param name="Syz">终点桩号</param>
        /// <param name="T">竖曲线长度</param>
        /// <param name="R">竖曲线半径</param>
        /// <param name="H">变坡点高程</param>
        /// <param name="i1">后坡坡度（弧度制）</param>
        /// <param name="i2">前坡坡度（弧度制）</param>
        /// <param name="sig">坡度的凹凸（1：凸，-1：凹）</param>
        /// <returns></returns>
        public static double HightCompute(double Kp, double K, double Szy, double Syz,
            double T, double R, double H, double i1, double i2, double sig)
        {
            double x;           //任一点桩号与竖曲线起点桩号之差
            double y, h1, h = 0, d;
            if (sig == 1 || sig == -1)
            {
                if (Kp < K && Kp >= Szy)
                {
                    x = Kp - Szy;
                    y = x * x / (2 * R);
                    d = Math.Abs(K - Kp);
                    h1 = H - d * i1;
                    h = h1 - sig * y;
                }
                else if (Kp < K && Kp < Szy)
                {
                    d = Math.Abs(K - Kp);
                    h1 = H - d * i1;
                    h = h1;
                }
                else if (Kp > K && Kp < Syz)
                {
                    x = Syz - Kp;
                    y = x * x / (2 * R);
                    d = Math.Abs(K - Kp);
                    h1 = H + d * i2;
                    h = h1 - sig * y;
                }
                else if (Kp > K && Kp > Syz)
                {
                    d = Math.Abs(K - Kp);
                    h1 = H + d * i2;
                    h = h1;
                }
            }
            else
                throw new Exception("Tools:HightCompute:请输入正确的参数");
            return h;
        }
        #endregion

        #region 坐标转换
        /// <summary>
        /// 大地坐标转空间直角坐标
        /// </summary>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        /// <param name="H">大地高</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="X">空间直角X坐标</param>
        /// <param name="Y">空间直角Y坐标</param>
        /// <param name="Z">空间直角Z坐标</param>
        public static void Geo2Comp(double B, double L, double H, double a, double e,
            out double X, out double Y, out double Z)
        {
            double dRb, dRl, N;
            //度分秒转弧度制
            dRb = Dms2Rad(B);
            dRl = Dms2Rad(L);
            N = a / Math.Sqrt(1 - e * e * Math.Sin(dRb) * Math.Sin(dRb));
            X = (N + H) * Math.Cos(dRb) * Math.Cos(dRl);
            Y = (N + H) * Math.Cos(dRb) * Math.Sin(dRl);
            Z = (N * (1 - e * e) + H) * Math.Sin(dRb);
        }
        /// <summary>
        /// 空间直角坐标转换至大地坐标
        /// </summary>
        /// <param name="X">空间直角X坐标</param>
        /// <param name="Y">空间直角Y坐标</param>
        /// <param name="Z">空间直角Z坐标</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        /// <param name="H">大地高</param>
        public static void Comp2Geo(double X, double Y, double Z, double a, double e,
            out double B, out double L, out double H)
        {
            double h, pcs, N;
            const double eps = 1e-5;//设定精度
            L = (Math.Abs(X) >= eps) ? Math.Atan(Y / X) : 0;
            //判别点落在赤道平面上三四象限上情况并加以修正
            if (X < 0)
                L = Y >= 0 ? (L + Math.PI) : (L - Math.PI);
            if ((Math.Abs(Y) >= eps) || (Math.Abs(X) >= eps))  //计算纬度和高程
            {
                B = Math.Atan(Z / Math.Sqrt(X * X + Y * Y));    //设定纬度初始值进行迭代
                H = 0;
                do
                {
                    N = a / Math.Sqrt(1 - e * e * Math.Sin(B) * Math.Sin(B));
                    h = Math.Sqrt(X * X + Y * Y) / Math.Cos(B) - N;
                    B = Math.Atan(Z / (Math.Sqrt(X * X + Y * Y) * (1 + e * e * N / (N + h))));
                    pcs = h - H;
                    H = h;
                } while (Math.Abs(pcs) >= eps);
            }
            else                               //单独定义Z轴上的大地坐标
            {
                B = 0;
                N = a / Math.Sqrt(1 - e * e * Math.Sin(B) * Math.Sin(B));
                H = Math.Sqrt(X * X + Y * Y) / Math.Cos(B) - N;
            }
            //必要修正，防止出现“”
            B += 1e-10;
            L += 1e-10;
            //弧度转化成度分秒
            B = Rad2Dms(B);
            L = Rad2Dms(L);
        }
        /// <summary>
        /// 高斯正算
        /// </summary>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        /// <param name="L0">中央经线（度分秒）</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="x">高斯投影X坐标</param>
        /// <param name="y">高斯投影Y坐标</param>
        public static void GaussForward(double B, double L, double L0, double a, double e,
            out double x, out double y)
        {
            B = Dms2Rad(B);
            L = Dms2Rad(L);
            double b = a * Math.Sqrt(1 - e * e);
            double W = Math.Sqrt(1 - e * e * Math.Sin(B) * Math.Sin(B));
            double N = a / W;
            double M = a * (1 - e * e) / Math.Pow(W, 3);
            double t = Math.Tan(B);
            double eitef = (a * a - b * b) * Math.Cos(B) * Math.Cos(B) / (b * b);
            double l = Dms2Rad(L) - Dms2Rad(L0);
            //主曲率半径计算
            double m0, m2, m4, m6, m8, n0, n2, n4, n6, n8;
            m0 = a * (1 - e * e); n0 = a;
            m2 = 3.0 / 2.0 * e * e * m0; n2 = 1.0 / 2.0 * e * e * n0;
            m4 = 5.0 / 4.0 * e * e * m2; n4 = 3.0 / 4.0 * e * e * n2;
            m6 = 7.0 / 6.0 * e * e * m4; n6 = 5.0 / 6.0 * e * e * n4;
            m8 = 9.0 / 8.0 * e * e * m6; n8 = 7.0 / 8.0 * e * e * n6;
            //子午线曲率半径
            double a0, a2, a4, a6, a8;
            a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
            a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7.0 / 16.0 * m8;
            a4 = m4 / 8.0 + 3.0 / 16.0 * m6 + 7.0 / 32.0 * m8;
            a6 = m6 / 32 + m8 / 16;
            a8 = m8 / 128;

            double X = a0 * B - a2 / 2 * Math.Sin(2 * B) + a4 / 4 * Math.Sin(4 * B) - a6 / 6 * Math.Sin(6 * B) + a8 / 8 * Math.Sin(8 * B);
            x = X + N / 2 * t * Math.Cos(B) * Math.Cos(B) * l * l +
                N / 24 * t * (5 - t * t + 9 * eitef + 4 * Math.Pow(eitef, 2)) *
                Math.Pow(Math.Cos(B), 4) * Math.Pow(l, 4) + N / 720 * t *
                (61 - 58 * t * t + Math.Pow(t, 4)) * Math.Pow(Math.Cos(B), 6) * Math.Pow(l, 6);
            y = N * Math.Cos(B) * l + N / 6 * (1 - t * t + eitef) * Math.Pow(Math.Cos(B), 3) *
                Math.Pow(l, 3) + N / 120 *
                (5 - 18 * t * t + Math.Pow(t, 4) + 14 * eitef - 58 * eitef * t * t) *
                Math.Pow(Math.Cos(B), 5) * Math.Pow(l, 5);
        }
        /// <summary>
        /// 计算子午线弧长
        /// </summary>
        /// <param name="B">大地纬度（弧度）</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e1">椭球第一偏心率</param>
        /// <returns></returns>
        public static double ComptX(double B, double a, double e1)
        {
            int N = 5;          //计算项数
            double sb = Math.Cos(B);
            double sb2 = sb * sb;
            double scb = sb * Math.Cos(B);
            double e = 1.0;
            double In = B;
            double sum = B;
            double k = 0.0;
            for (int n = 1; n <= N; ++n)
            {
                k = 0.5 / n;
                e = e * (1.0 + k) * e1 * e1;
                In = In - k * (In + scb);
                sum = sum + e * In;
                scb = scb * sb2;
            }
            return a * (1.0 - e1 * e1) * sum;
        }
        /// <summary>
        /// 分带高斯正算
        /// </summary>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        /// <param name="Zone">度带（3：三度带）</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="dX">高斯投影X坐标</param>
        /// <param name="dY">高斯投影Y坐标</param>
        public static void GaussZoneForward(double B, double L, double Zone, double a, double e,
            out double dX, out double dY)
        {
            L = Dms2Deg(L);
            if (Zone == 3)      //判断度带情况
                L -= 1.5;
            int ZoneNum = (int)(L / Zone) + 1;
            if (ZoneNum <= 0)
                ZoneNum = (int)(360 / Zone) + ZoneNum;  //西半球情况
            double L0 = (ZoneNum - 0.5) * Zone + 1.5;  //计算中央经度
            GaussForward(B, L, Deg2Dms(L0), a, e, out dX, out dY); //高斯正算
            dY = (dY >= 0 ? dY + ZoneNum * 1e-6 : dY - ZoneNum * 1e-6);   //添加带号
        }
        /// <summary>
        /// 高斯反算
        /// </summary>
        /// <param name="x">高斯投影X坐标</param>
        /// <param name="y">高斯投影Y坐标</param>
        /// <param name="L0">中央经线（度分秒）</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        public static void GaussInverse(double x, double y, double L0, double a, double e,
            out double B, out double L)
        {
            double b = a * Math.Sqrt(1 - e * e);
            double m0, m2, m4, m6, m8, n0, n2, n4, n6, n8;
            m0 = a * (1 - e * e); n0 = a;
            m2 = 3.0 / 2.0 * e * e * m0; n2 = 1.0 / 2.0 * e * e * n0;
            m4 = 5.0 / 4.0 * e * e * m2; n4 = 3.0 / 4.0 * e * e * n2;
            m6 = 7.0 / 6.0 * e * e * m4; n6 = 5.0 / 6.0 * e * e * n4;
            m8 = 9.0 / 8.0 * e * e * m6; n8 = 7.0 / 8.0 * e * e * n6;

            //子午线曲率半径
            double a0, a2, a4, a6, a8;
            a0 = m0 + m2 / 2 + 3 * m4 / 8 + 5 * m6 / 16 + 35 * m8 / 128;
            a2 = m2 / 2 + m4 / 2 + 15 * m6 / 32 + 7.0 / 16.0 * m8;
            a4 = m4 / 8.0 + 3.0 / 16.0 * m6 + 7.0 / 32.0 * m8;
            a6 = m6 / 32 + m8 / 16;
            a8 = m8 / 128;

            double X = x;
            double FBf = 0;
            double Bf0 = X / a0, Bf1 = 0;
            while ((Bf0 - Bf1) >= 0.0001)
            {
                Bf1 = Bf0;
                FBf = a0 * Bf0 - a2 / 2 * Math.Sin(2 * Bf0) + a4 / 4 * Math.Sin(4 * Bf0) - a6 / 6 * Math.Sin(6 * Bf0) + a8 / 8 * Math.Sin(8 * Bf0);
                Bf0 = (X - FBf) / a0;
            }
            double Bf = Bf0;
            double Vf = Math.Sqrt(1 + (a * a - b * b) * Math.Cos(Bf) * Math.Cos(Bf) / (b * b));
            double tf = Math.Tan(Bf);
            double Nf = a / Math.Sqrt(1 - Math.Pow(e * Math.Sin(Bf), 2));
            double eiteffang = (a * a - b * b) * Math.Cos(Bf) * Math.Cos(Bf) / (b * b);
            double Bdu = Rad2Deg(Bf) - 1 / 2.0 * Vf * Vf * tf *
                (Math.Pow((y / Nf), 2) - 1.0 / 12 * (5 + 3 * tf * tf + eiteffang - 9 * eiteffang * tf * tf) *
                Math.Pow((y / Nf), 4) + 1.0 / 360.0 * (61 + 90 * tf * tf + 45 * tf * tf) *
                Math.Pow((y / Nf), 6)) * 180 / Math.PI;
            double ldu = 1.0 / Math.Cos(Bf) * (y / Nf + 1.0 / 6.0 * (1 + 2 * tf * tf + eiteffang) *
                Math.Pow((y / Nf), 3) + 1.0 / 120.0 *
                (5 + 28 * tf * tf + 24 * tf * tf + 6 * eiteffang + 8 * eiteffang * tf * tf) *
                Math.Pow((y / Nf), 5)) * 180.0 / Math.PI;

            B = Deg2Dms(Bdu);
            L = Rad2Dms(Dms2Rad(L0) + Deg2Rad(ldu));
        }
        /// <summary>
        /// 分带高斯反算
        /// </summary>
        /// <param name="x">高斯投影X坐标</param>
        /// <param name="y">高斯投影Y坐标</param>
        /// <param name="Zone">度带（3：3度带）</param>
        /// <param name="a">椭球长半轴</param>
        /// <param name="e">椭球偏心率 第一偏心率</param>
        /// <param name="B">大地纬度（度分秒）</param>
        /// <param name="L">大地经度（度分秒）</param>
        public static void GaussZoneInverse(double x, double y, double Zone, double a, double e,
            out double B, out double L)
        {
            int ZoneNum = (int)(y / 1e6);           //提取带号
            y = y - ZoneNum * 1e6;              //计算去带号和改正值后的Y左边
            if (ZoneNum > 180 / Zone)           //西半球情况
                ZoneNum = ZoneNum - (int)(360 / Zone);
            double L0 = (ZoneNum - 0.5) * Zone;
            //判断度带情况
            if (Zone == 3)
                L0 += 1.5;
            GaussInverse(x, y, Deg2Dms(L0), a, e, out B, out L);
        }
        /// <summary>
        /// 平面点
        /// </summary>
        public struct HPOINT
        {
            /// <summary>
            /// 该点的X坐标
            /// </summary>
            public double X;
            /// <summary>
            /// 该点的Y坐标
            /// </summary>
            public double Y;
        }
        /// <summary>
        /// 三维点
        /// </summary>
        public struct SPOINT
        {
            /// <summary>
            /// 该点的X坐标
            /// </summary>
            public double X;
            /// <summary>
            /// 该点的Y坐标
            /// </summary>
            public double Y;
            /// <summary>
            /// 该点的Z坐标
            /// </summary>
            public double Z;
        }
        /// <summary>
        /// 转换四参数
        /// </summary>
        public struct FourPara
        {
            /// <summary>
            /// X轴变化值
            /// </summary>
            public double DeltaX;
            /// <summary>
            /// Y轴变化值
            /// </summary>
            public double DeltaY;
            /// <summary>
            /// 尺度参数
            /// </summary>
            public double Scale;
            /// <summary>
            /// 旋转参数
            /// </summary>
            public double Sita;
        }
        /// <summary>
        /// 转换七参数
        /// </summary>
        public struct SevenPara
        {
            /// <summary>
            /// X轴旋转量（弧度）
            /// </summary>
            public double Alpha;
            /// <summary>
            ///  Y轴旋转量（弧度）
            /// </summary>
            public double Betta;
            /// <summary>
            ///  Z轴旋转量（弧度）
            /// </summary>
            public double Gamma;
            /// <summary>
            /// X轴变化值
            /// </summary>
            public double DeltaX;
            /// <summary>
            /// Y轴变化值
            /// </summary>
            public double DeltaY;
            /// <summary>
            /// Z轴变化值
            /// </summary>
            public double DeltaZ;
            /// <summary>
            /// 尺度参数
            /// </summary>
            public double Scale;
        }
        /// <summary>
        /// 两点法求平面转换四参数
        /// </summary>
        /// <param name="PointS">转换前两点坐标</param>
        /// <param name="PointT">转换后两点坐标</param>
        /// <returns></returns>
        public static FourPara ComFourParaUsingTwoPoints(HPOINT[] PointS, HPOINT[] PointT)
        {
            FourPara FP = new FourPara();
            FP.Sita = Azimuth(PointT[0].X, PointT[0].Y, PointT[1].X, PointT[1].Y) -
                Azimuth(PointS[0].X, PointS[0].Y, PointS[1].X, PointS[1].Y);
            FP.Scale = Dist(PointT[0].X, PointT[0].Y, PointT[1].X, PointT[1].Y) /
                Dist(PointS[0].X, PointS[0].Y, PointS[1].X, PointS[1].Y);
            FP.DeltaX = PointT[0].X - FP.Scale * Math.Cos(FP.Sita) * PointS[0].X +
                FP.Scale * Math.Sin(FP.Sita) * PointS[0].Y;
            FP.DeltaY = PointT[0].Y - FP.Scale * Math.Sin(FP.Sita) * PointS[0].X -
                FP.Scale * Math.Cos(FP.Sita) * PointS[0].Y;
            return FP;
        }
        /// <summary>
        /// 多点法求平面转换四参数
        /// </summary>
        /// <param name="PointS">转换前的点</param>
        /// <param name="PointT">转换后的点</param>
        /// <param name="PointCount">点数</param>
        /// <returns></returns>
        public static FourPara ComFourPara(HPOINT[] PointS, HPOINT[] PointT, int PointCount)
        {
            FourPara FP = new FourPara();
            double u, v, K, Sita, DeltaX, DeltaY;
            //初始化
            DeltaX = 0;
            DeltaY = 0;
            u = 1;
            v = 0;
            int intCount = PointCount;
            CMatrix dx = new CMatrix(4, 1);                 //待估参数改正数
            CMatrix B = new CMatrix(2 * intCount, 4);       //误差方程系数矩阵
            CMatrix W = new CMatrix(2 * intCount, 1);       //误差方程常数项
            CMatrix BT, N, InvN, BTW;

            for (int i = 0; i < intCount; i++)             //计算误差方程系数矩阵
            {
                B[2 * i, 0] = 1;
                B[2 * i, 1] = 0;
                B[2 * i, 2] = PointS[i].X;
                B[2 * i, 3] = -PointS[i].Y;

                B[2 * i + 1, 0] = 0;
                B[2 * i + 1, 1] = 1;
                B[2 * i + 1, 2] = PointS[i].Y;
                B[2 * i + 1, 3] = PointS[i].X;
            }
            for (int i = 0; i < intCount; i++)             //计算误差方程系常数
            {
                W[2 * i, 0] = PointT[i].X - u * PointS[i].X + v * PointS[i].Y - DeltaX;
                W[2 * i + 1, 0] = PointT[i].Y - u * PointS[i].Y - v * PointS[i].X - DeltaY;
            }
            //最小二乘求解
            BT = B.Translocation();
            N = BT * B;
            InvN = N.Inv();
            BTW = BT * W;
            dx = InvN * BTW;
            DeltaX = DeltaX + dx[0, 0];
            DeltaY = DeltaY + dx[1, 0];
            u = u + dx[2, 0];
            v = v + dx[3, 0];

            Sita = Math.Atan(v / u);
            K = u / Math.Cos(Sita);
            FP.DeltaX = DeltaX;
            FP.DeltaY = DeltaY;
            FP.Scale = K;
            FP.Sita = Sita;
            return FP;
        }
        /// <summary>
        /// 多点法求平面转换七参数
        /// </summary>
        /// <param name="PointS">转换前的点</param>
        /// <param name="PointT">转换后的点</param>
        /// <param name="PointCount">点数</param>
        /// <returns></returns>
        public static SevenPara ComSevenPara(SPOINT[] PointS, SPOINT[] PointT, int PointCount)
        {
            SevenPara SP = new SevenPara();
            CMatrix B = new CMatrix(PointCount * 3, 7); //如果是个已知点，9*7矩阵V=B*X-L中的矩阵B
            CMatrix dX = new CMatrix(7, 1);             //V=B*X-L中的矩阵X
            CMatrix L = new CMatrix(PointCount * 3, 1); //V=B*X-L中的矩阵L，如果是个已知点，9*1矩阵
            CMatrix BT, N, InvN, BTL;

            for (int i = 0; i < PointCount * 3; i++)
                if (i % 3 == 0)
                    L[i, 0] = PointT[i / 3].X;
                else if (i % 3 == 1)
                    L[i, 0] = PointT[i / 3].Y;
                else if (i % 3 == 2)
                    L[i, 0] = PointT[i / 3].Z;
            //B矩阵
            for (int i = 0; i < PointCount * 3; i++)
                if (i % 3 == 0)
                {
                    B[i, 0] = 1;
                    B[i, 1] = 0;
                    B[i, 2] = 0;
                    B[i, 3] = PointS[i / 3].X;
                    B[i, 4] = 0;
                    B[i, 5] = -PointS[i / 3].Z;
                    B[i, 6] = PointS[i / 3].Y;
                }
                else if (i % 3 == 1)
                {
                    B[i, 0] = 0;
                    B[i, 1] = 1;
                    B[i, 2] = 0;
                    B[i, 3] = PointS[i / 3].Y;
                    B[i, 4] = PointS[i / 3].Z;
                    B[i, 5] = 0;
                    B[i, 6] = -PointS[i / 3].X;
                }
                else if (i % 3 == 2)
                {
                    B[i, 0] = 0;
                    B[i, 1] = 0;
                    B[i, 2] = 1;
                    B[i, 3] = PointS[i / 3].Z;
                    B[i, 4] = -PointS[i / 3].Y;
                    B[i, 5] = PointS[i / 3].X;
                    B[i, 6] = 0;
                }
            BT = B.Translocation();
            N = BT * B;
            InvN = N.Inv();
            BTL = BT * L;
            dX = InvN * BTL;
            //结果
            SP.DeltaX = dX[0, 0];
            SP.DeltaY = dX[1, 0];
            SP.DeltaZ = dX[2, 0];
            SP.Scale = dX[3, 0] - 1;
            SP.Alpha = dX[4, 0] / dX[3, 0];
            SP.Betta = dX[5, 0] / dX[3, 0];
            SP.Gamma = dX[6, 0] / dX[3, 0];
            return SP;
        }
        /// <summary>
        /// 七参数空间直角坐标转换
        /// </summary>
        /// <param name="PointSource">转换前坐标</param>
        /// <param name="SP">转换参数</param>
        /// <returns></returns>
        public static SPOINT SPointsTrans(SPOINT PointSource, SevenPara SP)
        {
            SPOINT PointTarget = new SPOINT();
            PointTarget.X = (1 + SP.Scale) *
                (PointSource.X + SP.Gamma * PointSource.Y - SP.Betta * PointSource.Z) + SP.DeltaX;
            PointTarget.Y = (1 + SP.Scale) *
                (-SP.Gamma * PointSource.X + PointSource.Y + SP.Alpha * PointSource.Z) + SP.DeltaY;
            PointTarget.Z = (1 + SP.Scale) *
                (SP.Betta * PointSource.X - SP.Alpha * PointSource.Y + PointSource.Z) + SP.DeltaZ;
            return PointTarget;
        }
        /// <summary>
        /// 四参数平面直角坐标转换
        /// </summary>
        /// <param name="PointSource">转换前坐标</param>
        /// <param name="FP">转换参数</param>
        /// <returns></returns>
        public static HPOINT HPointTrans(HPOINT PointSource, FourPara FP)
        {
            HPOINT PointTarget = new HPOINT();
            PointTarget.X = FP.DeltaX + FP.Scale * PointSource.X * Math.Cos(FP.Sita) -
                FP.Scale * PointSource.Y * Math.Sin(FP.Sita);
            PointTarget.Y = FP.DeltaY + FP.Scale * PointSource.X * Math.Sin(FP.Sita) +
                FP.Scale * PointSource.Y * Math.Cos(FP.Sita);
            return PointTarget;
        }
        #endregion

    }
}
