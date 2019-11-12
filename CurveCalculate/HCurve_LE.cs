using System;

namespace EngineeringSurveyTools.CurveCalculate
{
    /// <summary>
    /// 平曲线线元
    /// </summary>
    class HCurve_LE
    {
        /// <summary>
        /// 线型（直线：0，直缓圆：1，圆缓直：2，圆弧：3）
        /// </summary>
        public int type;
        /// <summary>
        /// 起始里程
        /// </summary>
        public double K0;
        /// <summary>
        /// 结束里程
        /// </summary>
        public double K1;
        /// <summary>
        /// 起始方位角
        /// </summary>
        public double azimuth;
        /// <summary>
        /// 起始北坐标
        /// </summary>
        public double x0;
        /// <summary>
        /// 起始东坐标
        /// </summary>
        public double y0;
        /// <summary>
        /// 转弯半径
        /// </summary>
        public double R;
        /// <summary>
        /// 曲线长
        /// </summary>
        public double L;
        /// <summary>
        /// 起点曲率
        /// </summary>
        public double rhoA;
        /// <summary>
        /// 终点曲率
        /// </summary>
        public double rhoB;
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="type">线型（直线：0，直缓圆：1，圆缓直：2，圆弧：3）</param>
        /// <param name="K0">起始里程</param>
        /// <param name="A">起始方位角</param>
        /// <param name="x0">起始北坐标</param>
        /// <param name="y0">起始东坐标</param>
        /// <param name="R">转弯半径</param>
        /// <param name="L">平曲线线元长</param>
        public HCurve_LE(int type,double K0,double A,double x0,double y0,double R,double L)
        {
            this.type = type;
            this.K0 = K0;
            azimuth = A;
            this.x0 = x0;
            this.y0 = y0;
            this.R = R;
            this.L = L;
            K1 = K0 + L;
            if (R == 0)
                this.R = double.MaxValue;
            switch (type)
            {
                case 0:
                    {
                        rhoA = 0;
                        rhoB = 0;
                        break;
                    }
                case 1:
                    {
                        rhoA = 0;
                        rhoB = 1 / R;
                        break;
                    }
                case 2:
                    {
                        rhoA = 1 / R;
                        rhoB = 0;
                        break;
                    }
                case 3:
                    {
                        rhoA = 1 / R;
                        rhoB = 1 / R;
                        break;
                    }
                default: throw new Exception("平曲线线元线型未知！");
            }
        }
    }
}
