using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringSurveyTools.CurveCalculate
{
    enum BrokenType
    {
        LongChainage = 0,   //长链
        ShortChainage = 1,  //短链
        NumSwitch = 2,      //冠号切换
        Spare1 = 3,         //空余位1
        Spare2 = 4          //空余位2
    };
    class BrokenChainage
    {
        public string name;             //点名
        public BrokenType type;         //点类型
        public double beforeK;          //断链前里程，非负数（实际里程）
        public double afterK;           //断链后里程，非负数（实际里程）
        public double L;                //断链长度
        public string prefix;        //断链前冠号
        public string suffix;         //断链后冠号
        public double continuousK;      //断链点的连续里程
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="type">断链点类型(0:长链，1:短链，2:冠号切换)</param>
        /// <param name="bK">断链前里程，非负值</param>
        /// <param name="aK">断链后里程，非负值</param>
        /// <param name="pre">断链前冠号</param>
        /// <param name="suf">断链后冠号</param>
        /// <param name="name">断链点名</param>
        public BrokenChainage(int type,double bK,double aK,
            string pre,string suf,string name="")
        {
            if (bK < 0 || aK < 0)
                throw new Exception("BrokenChainage:BrokenChainage:断链里程为负值");
            this.type = (BrokenType)type;
            beforeK = bK;
            afterK = aK;
            L = bK - aK;
            if (((L > 0) && (type != 0)) || ((L < 0) && (type != 1)) ||
                ((L == 0) && (type != 2)))
                throw new Exception("BrokenChainage:BrokenChainage:断链类型与长度值不匹配");
            prefix = pre;
            suffix = suf;
        }


    }
}
