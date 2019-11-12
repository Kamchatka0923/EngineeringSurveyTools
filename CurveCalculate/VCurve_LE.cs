using System;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 竖曲线线元
    /// </summary>
    class VCurve_LE
    {
        public double i1;           //后坡坡度
        public double i2;           //前坡坡度
        public double K;            //变坡点里程
        public double H;            //变坡点高程
        public double R;            //半径
        //竖曲线要素
        private double L;            //竖曲线长度
        private double T;            //切线长度
        private double E;            //竖曲线外距
        private double Szy;          //ZY点里程
        private double Syz;          //YZ点里程

        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="K">里程</param>
        /// <param name="H">高程</param>
        /// <param name="R">半径</param>
        public VCurve_LE(double K, double H, double R)
        {
            this.K = K;
            this.H = H;
            this.R = R;
            i1 = 999;
            i2 = 999;
        }
        /// <summary>
        /// 添加坡度信息
        /// </summary>
        /// <param name="i1">后坡坡度</param>
        /// <param name="i2">前坡坡度</param>
        public void Add_i(double i1, double i2)
        {
            this.i1 = i1;
            this.i2 = i2;
            if (i1 == i2)
                return;
            ParaOfVerticurve(i1, i2, R, K, H, out L, out T, out E, out Szy, out Syz);
        }
        /// <summary>
        /// 根据里程计算高程
        /// </summary>
        /// <param name="Kp">未知高程点里程</param>
        /// <returns></returns>
        public double HightComputeByK(double Kp)
        {
            if (i1 == 999 || i2 == 999)
                throw new Exception("VerticalCurve:HightCompute:未添加坡度信息");
            if (i1 == i2)
                throw new Exception("VerticalCurve:HightCompute:该竖曲线要素无法计算");
            double sig = 1;
            if ((i1 - i2) < 0)
                sig = -1;
            double h = HightCompute(Kp, K, Szy, Syz, T, R, H, i1, i2, sig);
            return h;
        }
        /// <summary>
        /// 计算竖曲线要素
        /// </summary>
        /// <param name="i1">后坡坡度</param>
        /// <param name="i2">前坡坡度</param>
        /// <param name="R">竖曲线半径</param>
        /// <param name="K">变坡点桩号</param>
        /// <param name="H">变坡点高程</param>
        /// <param name="L">竖曲线长度</param>
        /// <param name="T">切线长度</param>
        /// <param name="E">竖曲线外距</param>
        /// <param name="Szy">起点桩号</param>
        /// <param name="Syz">终点桩号</param>
        private void ParaOfVerticurve(double i1, double i2, double R, double K, double H,
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
        private double HightCompute(double Kp, double K, double Szy, double Syz,
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
    }
}
