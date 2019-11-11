using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SurveyTools
{
    /// <summary>
    /// 矩阵类
    /// </summary>
    public class CMatrix
    {
        private int iRow;
        private int iCol;
        private double[][] dMatData;
        /// <summary>
        /// 矩阵的行数
        /// </summary>
        /// <returns></returns>
        public int Row() { return iRow; }
        /// <summary>
        /// 矩阵的列数
        /// </summary>
        /// <returns></returns>
        public int Col() { return iCol; }
        /// <summary>
        /// 拷贝构造函数
        /// </summary>
        /// <param name="m"></param>
        public CMatrix(CMatrix m)
        {
            iRow = m.Row();
            iCol = m.Col();
            dMatData = new double[iRow][];
            for (int i = 0; i < iRow; i++)
            {
                dMatData[i] = new double[iCol];
                for (int j = 0; j < iCol; j++)
                    Array.Copy(m.dMatData[i], dMatData[i], iCol);
            }
        }
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="row">矩阵的行数</param>
        /// <param name="col">矩阵的列数</param>
        public CMatrix(int row, int col)
        {
            iRow = row;
            iCol = col;
            dMatData = new double[row][];
            for (int i = 0; i < row; i++)
                dMatData[i] = new double[col];
        }
        /// <summary>
        /// 构造函数（矩阵行列数未知）
        /// </summary>
        public CMatrix()
        {
            dMatData = new double[0][];
        }
        /// <summary>
        /// 矩阵索引
        /// </summary>
        /// <param name="row">行</param>
        /// <param name="col">列</param>
        /// <returns></returns>
        public double this[int row, int col]
        {
            get
            {
                if (row >= iRow || col <= iCol)
                    throw new Exception("CMatrix:[,]:Index out of range!");
                return dMatData[row][col];
            }
            set
            {
                dMatData[row][col] = value;
            }
        }
        /// <summary>
        /// 重载+
        /// </summary>
        public static CMatrix operator +(CMatrix m1, CMatrix m2)
        {
            if ((m1.Col() != m2.Col()) || (m1.Row() != m2.Row()))
                throw new Exception
                    ("Cmatrix:operator+:The two matrix have difference size!");
            CMatrix matTmp = new CMatrix(m1.Row(), m1.Col());
            for (int i = 0; i < m1.Row(); i++)
                for (int j = 0; j < m1.Col(); j++)
                    matTmp[i, j] = m1[i, j] + m2[i, j];
            return matTmp;
        }
        /// <summary>
        /// 重载-
        /// </summary>
        public static CMatrix operator -(CMatrix m1, CMatrix m2)
        {
            if ((m1.Col() != m2.Col()) || (m1.Row() != m2.Row()))
                throw new Exception
                    ("Cmatrix:operator-:The two matrix have difference size!");
            CMatrix matTmp = new CMatrix(m1.Row(), m1.Col());
            for (int i = 0; i < m1.Row(); i++)
                for (int j = 0; j < m1.Col(); j++)
                    matTmp[i, j] = m1[i, j] - m2[i, j];
            return matTmp;
        }
        /// <summary>
        /// 重载*
        /// </summary>
        public static CMatrix operator *(CMatrix m1, CMatrix m2)
        {
            if ((m1.Col() != m2.Row()) || (m1.Row() != m2.Col()))
                throw new Exception
                    ("Cmatrix:operator*:The two matrix have difference size!");
            CMatrix matTmp = new CMatrix(m1.Row(), m2.Col());
            for (int i = 0; i < m1.Row(); i++)
                for (int j = 0; j < m2.Col(); j++)
                {
                    double sum = 0;
                    for (int k = 0; k < m2.Row(); k++)
                        sum += m1[i, k] * m2[k, j];
                    matTmp[i, j] = sum;
                }
            return matTmp;
        }
        /// <summary>
        /// 矩阵转置
        /// </summary>
        /// <returns></returns>
        public CMatrix Translocation()
        {
            CMatrix matTmp = new CMatrix(this.Col(), this.Row());
            for (int i = 0; i < this.Row(); i++)
                for (int j = 0; j < this.Col(); j++)
                    matTmp[j, i] = this[i, j];
            return matTmp;
        }

        /// <summary>
        /// 调整矩阵大小，原有值不变
        /// </summary>
        public void SetSize(int row, int col)
        {
            if (row == iRow && col == iCol)
                return;
            double[][] rsData = new double[row][];
            for (int i = 0; i < row; i++)
            {
                rsData[i] = new double[col];
                for (int j = 0; j < col; j++)
                    rsData[i][j] = 0;
            }
            int minRow = (iRow > row) ? row : iRow;
            int minCol = (iCol > col) ? col : iCol;
            for (int i = 0; i < minRow; i++)
                Array.Copy(dMatData[i], rsData[i], minCol);
            dMatData = rsData;
            iRow = row;
            iCol = col;
            return;
        }
        /// <summary>
        /// 矩阵求逆
        /// </summary>
        /// <returns></returns>
        public CMatrix Inv()
        {
            if (iRow != iCol)
                throw new Exception("CMatrix:Inv:待求逆的矩阵行列不相等！");
            int i, j, k, vv;
            CMatrix InvMat = new CMatrix(iRow, iRow);
            InvMat = new CMatrix(this);
            int[] MainRow = new int[iRow];
            int[] MainCol = new int[iRow];  //用于记录主元素的行和列
            double dMainCell;               //主元元素的值
            double dTemp;                   //临时变量
            for (k = 0; k < iRow; k++)
            {
                dMainCell = 0;
                //全选主元
                for (i = k; i < iRow; i++)
                    for (j = k; j < iRow; j++)
                    {
                        dTemp = Math.Abs(InvMat[i, j]);
                        if (dTemp > dMainCell)
                        {
                            dMainCell = dTemp;
                            MainRow[k] = i;
                            MainCol[k] = j;
                        }
                    }

                if (Math.Abs(dMainCell) < 0.0000000000001)
                    throw new Exception("CMatrix:Inv:矩阵秩亏");//矩阵秩亏，不能求逆
                if (MainRow[k] != k)
                    for (j = 0; j < iRow; j++)      //交换行
                    {
                        vv = MainRow[k];
                        dTemp = InvMat[k, j];
                        InvMat[k, j] = InvMat[vv, j];
                        InvMat[vv, j] = dTemp;
                    }
                if (MainCol[k] != k)                //交换列
                    for (i = 0; i < iRow; i++)
                    {
                        vv = MainRow[k];
                        dTemp = InvMat[i, k];
                        InvMat[i, k] = InvMat[i, vv];
                        InvMat[i, vv] = dTemp;
                    }
                InvMat[k, k] = 1.0 / InvMat[k, k];  //计算乘数
                for (j = 0; j < iRow; j++)                 //计算主行
                    if (j != k)
                        InvMat[k, j] = InvMat[k, j] * InvMat[k, k];
                for (i = 0; i < iRow; i++)                 //消元
                    if (i != k)
                        for (j = 0; j < iRow; j++)
                            if (j != k)
                                InvMat[i, j] -= InvMat[i, k] * InvMat[k, j];
                for (i = 0; i < iRow; i++)
                    if (i != k)
                        InvMat[i, k] = -InvMat[i, k] * InvMat[k, k];
            }
            for (k = iRow - 1; k >= 0; k--)
            {
                if (MainCol[k] != k)       //交换行
                    for (j = 0; j < iRow; j++)
                    {
                        vv = MainCol[k];
                        dTemp = InvMat[k, j];
                        InvMat[k, j] = InvMat[vv, j];
                        InvMat[vv, j] = dTemp;
                    }
                if (MainRow[k] != k)       //交换行
                    for (i = 0; i < iRow; i++)
                    {
                        vv = MainCol[k];
                        dTemp = InvMat[i, k];
                        InvMat[i, k] = InvMat[i, vv];
                        InvMat[i, vv] = dTemp;
                    }
            }
            return InvMat;
        }
    }
}
