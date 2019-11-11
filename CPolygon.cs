using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SurveyTools
{
    /// <summary>
    /// 多边形类
    /// </summary>
    public class CPolygon
    {
        /// <summary>
        /// 平面点
        /// </summary>
        public struct HPOINT
        {
            /// <summary>
            /// 平面点的坐标
            /// </summary>
            public double x, y;
        }
        private int iPointCount;//多边形顶点个数
        private HPOINT[] pPointData;//保存多边形顶点数据的数组
        /// <summary>
        /// 构造函数，多边形的顶点个数为0
        /// </summary>
        public CPolygon()
        {
            pPointData = null;
            iPointCount = 0;
        }
        /// <param name="n">多边形顶点个数</param>
        public CPolygon(int n)
        {
            pPointData = new HPOINT[n];
            iPointCount = n;
        }
        /// <summary>
        /// 设置第pos顶点（从开始）数据
        /// </summary>
        /// <param name="pos">顶点序数</param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public bool SetPoint(int pos, double x, double y)
        {
            if (pos >= iPointCount)
                return false;
            pPointData[pos].x = x;
            pPointData[pos].y = y;
            return true;
        }
        /// <summary>
        /// 获取顶点坐标
        /// </summary>
        /// <param name="pos">顶点序数</param>
        /// <param name="p">目标顶点</param>
        /// <returns></returns>
        public bool GetPoint(int pos, ref HPOINT p)
        {
            if (pos >= iPointCount)
                return false;
            p = pPointData[pos];
            return true;
        }
        /// <summary>
        /// 设置多边形的大小（改变顶点个数）
        /// 如果新设置的顶点个数小于原来的顶点个数，则只保留前面的顶点个数
        /// 反之扩充，原有顶点数据不变
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public bool SetSize(int n)
        {
            if (n < 3)
                return false;
            HPOINT[] pTemp = new HPOINT[n];
            int min = (iPointCount < n) ? iPointCount : n;
            Array.Copy(pPointData, pTemp, min);
            pPointData = pTemp;
            iPointCount = n;
            return true;
        }
        /// <summary>
        /// 计算多边形面积
        /// </summary>
        public double Area()
        {
            double dArea = 0;
            for (int i = 0; i < iPointCount; i++)
            {
                if (i == iPointCount - 1)
                    dArea += 0.5 * (pPointData[0].x + pPointData[i].x) *
                        (pPointData[0].y - pPointData[i].y);
                else
                    dArea += 0.5 * (pPointData[i + 1].x + pPointData[i].x) *
                        (pPointData[i + 1].y - pPointData[i].y);
            }
            return dArea;
        }
        /// <summary>
        /// 计算并返回多边形周长
        /// </summary>
        public double Perimeter()
        {
            double dDist;
            double p = 0;
            for (int i = 0; i < iPointCount; i++)
            {
                if (i == iPointCount - 1)
                    dDist = CommonTools.Dist(pPointData[i].x, pPointData[i].y,
                        pPointData[0].x, pPointData[0].y);
                else
                    dDist = CommonTools.Dist(pPointData[i].x, pPointData[i].y,
                        pPointData[i + 1].x, pPointData[i + 1].y);
                p = p + dDist;
            }
            return p;
        }
    }
}
