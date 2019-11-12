using System.Collections.Generic;

namespace EngineeringSurveyTools.CurveCalculate
{
    class Milleage
    {
        public List<BrokenChainage> bC;//断链brokenChainages
        public List<double[]> K;

        /// <summary>
        /// 构造函数
        /// </summary>
        public Milleage()
        {
            bC = new List<BrokenChainage>();
            K = new List<double[]>();
        }
        /// <summary>
        /// 添加平曲线实际里程
        /// </summary>
        /// <param name="K0">起始里程</param>
        /// <param name="K1">结束里程</param>
        /// <param name="L">曲线长</param>
        public void AddHa(double K0, double K1, double L)
        {
            double[] Ki = new double[4];
            Ki[0] = K0;
            Ki[1] = K1;
            if (K.Count == 0)
                Ki[2] = K0;
            else
                Ki[2] = K[K.Count - 1][3];
            Ki[3] = Ki[2] + L;
            K.Add(Ki);
        }
        /// <summary>
        /// 添加平曲线实际里程
        /// </summary>
        /// <param name="K0">起始里程</param>
        /// <param name="L">曲线长</param>
        public void AddHa(double K0,double L)
        {
            double[] Ki = new double[4];
            Ki[0] = K0;
            if(K.Count==0)
                Ki[2] = K0;
            else
            {
                K[K.Count - 1][1] = K0;
                Ki[2] = K[K.Count - 1][3];
            }
            Ki[1] = Ki[0] + L;
            Ki[3] = Ki[2] + L;
            K.Add(Ki);
        }
        /// <summary>
        /// 添加竖曲线实际里程
        /// </summary>
        /// <param name="K0">起始里程</param>
        public void AddVa(double K0)
        {
            double[] Ki = new double[4];
            Ki[0] = K0;
            Ki[1] = Ki[2] = Ki[3] = K0;
            if (K.Count != 0)
            {
                K[K.Count-1][1] = K0;
                K[K.Count - 1][3] = K0;
            }
            K.Add(Ki);
        }
        /// <summary>
        /// 添加平曲线断链
        /// </summary>
        /// <param name="B">断链</param>
        public void AddHB(BrokenChainage B)
        {
            int n = GetIndexByB(B);
            if (n == -1)  //在里程外不添加
                return;
            bC.Add(B);
            if ((int)B.type == 0 || (int)B.type == 1)   //长链或短链
            {
                double[] k1 = new double[4];
                double[] k2 = new double[4];
                k1[0] = K[n][0];
                k1[1] = B.beforeK;
                k1[2] = K[n][2];
                k1[3] = B.beforeK - K[n][0] + K[n][2];
                k2[0] = B.afterK;
                k2[1] = K[n][1];
                k2[2] = k1[3];
                k2[3] = K[n][3];
                K.RemoveAt(n);
                K.Insert(n, k2);
                K.Insert(n, k1);
            }
            else        //其他类型不改变里程
                return;
        }
        /// <summary>
        /// 添加竖曲线断链
        /// </summary>
        /// <param name="B">断链</param>
        public void AddVB(BrokenChainage B)
        {
            if (K[K.Count-1][0] == K[K.Count - 1][1])
                K.RemoveAt(K.Count - 1);   //移除曲线长为0的线元（通常为最后一项）
            int n = GetIndexByB(B);
            if (n == -1)  //在里程外不添加
                return;
            bC.Add(B);
            if ((int)B.type == 0 || (int)B.type == 1)   //长链或短链
            {
                double[] k1 = new double[4];
                double[] k2 = new double[4];
                k1[0] = K[n][0];
                k1[1] = B.beforeK;
                k1[2] = K[n][0];
                k1[3] = B.beforeK;
                k2[0] = B.afterK;
                k2[1] = K[n][1];
                k2[2] = B.afterK;
                k2[3] = K[n][1];
                K.RemoveAt(n);
                K.Insert(n, k2);
                K.Insert(n, k1);
            }
            else        //其他断链类型不改变里程
                return;
            ContinueFlush();    //更新连续里程
        }
        /// <summary>
        /// 连续里程值更新
        /// </summary>
        private void ContinueFlush()
        {
            if (K.Count == 0 || K.Count == 1)   //没有里程或只有一个里程时无需更新连续里程
                return;
            for(int i=1;i<K.Count;i++)
            {
                double dl = K[i][1] - K[i][0];
                K[i][2] = K[i - 1][3];
                K[i][3] = K[i - 1][3] + dl;
            }
        }
        /// <summary>
        /// 得到断链所在的里程序号
        /// </summary>
        /// <param name="B">断链</param>
        /// <returns></returns>
        private int GetIndexByB(BrokenChainage B)
        {
            for (int i = 0; i < K.Count; i++)
                if (B.beforeK > K[i][0] && B.beforeK < K[i][1] && B.afterK > K[i][0] && B.afterK < K[i][1])
                    return i;
            return -1;
        }
        /// <summary>
        /// 根据实际里程计算连续里程（0为不在里程区间内）
        /// </summary>
        /// <param name="aK">实际里程(负数为重叠里程中后一个)</param>
        /// <returns></returns>
        public double GetContinueK(double aK)
        {
            double k1 = 0, k2 = 0;
            bool isAfter = false;
            if(aK<0)
            {
                aK = - aK;
                isAfter = true;
            }
            for (int i = 0; i < K.Count; i++)
                if (aK >= K[i][0] && aK <= K[i][1])
                {
                    k1 = aK - K[i][0] + K[i][2];
                    break;
                }
            for (int i = K.Count - 1; i >= 0; i--)
                if (aK >= K[i][0] && aK <= K[i][1])
                {
                    k2 = aK - K[i][0] + K[i][2];
                    break;
                }
            if (k1 == 0 && k2 == 0)
                return 0;
            if (isAfter)
                return k2;
            else
                return k1;
        }
        /// <summary>
        /// 根据连续里程计算实际里程(负数为重叠里程中后一个,0为不在里程区间内)
        /// </summary>
        /// <param name="cK">连续里程</param>
        /// <returns></returns>
        public double GetActualK(double cK)
        {
            double k1 = 0;
            for (int i = 0; i < K.Count; i++)
                if (cK > K[i][2] && cK <= K[i][3])
                {
                    k1 = cK - K[i][2] + K[i][0];
                    if (i > 0 && i < K.Count - 1)
                        if (k1 > K[i - 1][0] && k1 <= K[i - 1][1])
                            k1 = -k1;
                    break;
                }
            return k1;
        }
    }
}
