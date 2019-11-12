﻿using System;

namespace EngineeringSurveyTools.BothCorrect
{
    /// <summary>
    /// GB2007两化改正
    /// </summary>
    public class BothCorrectGB2007
    {
        private double DH;//归算到测区平均高程面上的测距边长度(m)
        private double DP;//测线的水平距离(m)
        private double HP;//测区的平均高程(m)
        private double Hm;//测距边两端点的平均高程(m)
        private double RA;//参考椭球体在测距边方向法截弧的曲率半径(m)
        private double D0;//归算到参考椭球面上的测距边长度(m)
        private double hm;//测区大地水准面高出参考椭球面的高差(m)
        private double Dg;//测距边在高斯投影面上的长度(m)
        private double ym;//测距边两端点横坐标的平均值(m)
        private double Rm;//测距边中点处在参考椭球面上的平均曲率半径(m)
        private double dy;//测距边两端点横坐标的增量(m)

        private double a;//参考椭球长半轴(m)
        private double b;//参考椭球短半轴(m)
        private double e;//参考椭球第一偏心率(m)
        private double e1;//参考椭球第二偏心率(m)
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="HP">测区的平均高程(m)</param>
        /// <param name="hm">测区大地水准面高出参考椭球面的高差(m)</param>
        /// <param name="a">参考椭球长半轴(m)</param>
        /// <param name="b">参考椭球短半轴(m)</param>
        public BothCorrectGB2007(double HP, double hm, double a, double b)
        {
            this.HP = HP;
            this.hm = hm;
            this.a = a;
            this.b = b;
            e = Math.Sqrt(a * a - b * b) / a;
            e1 = Math.Sqrt(a * a - b * b) / b;
        }
        /// <summary>
        /// 设置参考椭球体在测距边方向法截弧的曲率半径和测距边中点处在参考椭球面上的平均曲率半径
        /// </summary>
        /// <param name="B">大地纬度(弧度)</param>
        /// <param name="A">方位角(弧度)</param>
        private void SetRARm(double B, double A)
        {
            double M, N;
            M = a * (1 - e * e) * Math.Pow(1 - Math.Pow(e * Math.Sin(B), 2), -1.5);//子午圈曲率半径
            N = a * Math.Pow(1 - Math.Pow(e * Math.Sin(B), 2), -0.5);    //卯酉圈曲率半径
            RA = N / (1 + Math.Pow(e1 * Math.Cos(B) * Math.Cos(A), 2));
            Rm = Math.Sqrt(M * N);
        }
        /// <summary>
        /// 两化改正
        /// </summary>
        /// <param name="SA">A点到B点测线的水平距离(m)</param>
        /// <param name="SB">B点到A点测线的水平距离(m)</param>
        /// <param name="HA">A点高程(m)</param>
        /// <param name="HB">B点高程(m)</param>
        /// <param name="yA">A点横坐标(m)</param>
        /// <param name="yB">B点横坐标(m)</param>
        /// <param name="B">大地纬度(弧度)</param>
        /// <param name="A">方位角(弧度)</param>
        public double Calculate(double SA, double SB, double HA, double HB, double yA, double yB, double B, double A)
        {
            DP = (SA + SB) / 2;
            Hm = (HA + HB) / 2;
            ym = (yA + yB) / 2;
            dy = yB - yA;
            SetRARm(B, A);

            DH = DP * (1 + (HP - Hm) / RA);
            D0 = DP * (1 - (Hm + hm) / (RA + Hm + hm));
            Dg = D0 * (1 + ym * ym / (2 * Rm * Rm) + dy * dy / (24 * Rm * Rm));
            return Dg;
        }

    }
}
