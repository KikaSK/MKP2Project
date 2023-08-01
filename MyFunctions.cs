using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using OpenTK;

namespace MKP2___Template
{
    static class MyFunctions
    {
        private static List<int> GetIndices(int eU, int eV)
        {
            List<int> IndList = new List<int>();

            // indices for rectangles - quadruples  
            if (eU <= 0)
            {
                for (int i = 0; i < eV; i++)
                {
                    IndList.Add(i);
                    IndList.Add(i);
                    IndList.Add(i + 1);
                    IndList.Add(i + 1);
                }
            }
            else if (eV <= 0)
            {
                for (int i = 0; i < eU; i++)
                {
                    IndList.Add(i);
                    IndList.Add(i + 1);
                    IndList.Add(i + 1);
                    IndList.Add(i);
                }
            }
            else
            {
                for (int i = 0; i < eV; i++)
                    for (int j = 0; j < eU; j++)
                    {
                        IndList.Add(i * (eU + 1) + j);
                        IndList.Add(i * (eU + 1) + j + 1);
                        IndList.Add((i + 1) * (eU + 1) + j + 1);
                        IndList.Add((i + 1) * (eU + 1) + j);
                    }
            }
            return IndList;
        }
        public static Vector3 BilinearDeCasteljau(int m, int n, float u, float v, List<Vector3> ControlPoints)
        {
            List<int> Indices = GetIndices(m - 1, n - 1);
            if (ControlPoints.Count == 1) return ControlPoints[0];
            List<Vector3> newControlPoints = new List<Vector3>();
            for (int i = 0; i < Indices.Count / 4; ++i)
            {
                int i1 = Indices[4 * i];
                int i2 = Indices[4 * i + 1];
                int i3 = Indices[4 * i + 2];
                int i4 = Indices[4 * i + 3];
                Vector3 newPoint = ControlPoints[i1] * (1 - v) * (1 - u) + ControlPoints[i4] * (1 - u) * v + ControlPoints[i2] * (1 - v) * u + ControlPoints[i3] * u * v;
                newControlPoints.Add(newPoint);
            }
            return BilinearDeCasteljau(m - 1, n - 1, u, v, newControlPoints);
        }
        public static Vector3 GetUnitNormal(int m, int n, float u, float v, List<Vector3> ControlPoints)
        {
            List<Vector3> DifferenceU = new List<Vector3>();
            for (int i = 0; i < n - 1; ++i)
                for (int j = 0; j < m; ++j)
                    DifferenceU.Add(ControlPoints[(i + 1) * m + j] - ControlPoints[i * m + j]);

            List<Vector3> DifferenceV = new List<Vector3>();
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < m - 1; ++j)
                    DifferenceV.Add(ControlPoints[i * m + j + 1] - ControlPoints[i * m + j]);

            Vector3 dU = BilinearDeCasteljau(m, n - 1, u, v, DifferenceU);
            Vector3 dV = BilinearDeCasteljau(m - 1, n, u, v, DifferenceV);

            return Vector3.Cross(dU, dV).Normalized();
        }
        public static void GetSamplingAndNormals(ref List<Vector3> Sampling, ref List<Vector3> ControlPoints, int m, int n)
        {
            for (int i = 0; i < Sampling.Count; ++i)
            {
                var sample = Sampling[i];
                float u = (sample.X + 1) / 2;
                float v = (sample.Y + 1) / 2;
                Vector3 res = BilinearDeCasteljau(m, n, u, v, ControlPoints);
                Sampling[i] = res;
                //Vector3 normal = GetUnitNormal(m, n, u, v, ControlPoints);
                //Normals.Add(normal);
            }
        }
        public static Tuple<double, double, double> GetBC(Vector3 P, Vector3 A, Vector3 B, Vector3 C)
        {
            Point _P = new Point(P.X, P.Y);
            Point _A = new Point(A.X, A.Y);
            Point _B = new Point(B.X, B.Y);
            Point _C = new Point(C.X, C.Y);
            Vector AB = _B - _A;
            Vector AC = _C - _A;
            Vector PA = _A - _P;
            Vector PB = _B - _P;
            Vector PC = _C - _P;

            double areaABC = Vector.CrossProduct(AB, AC);
            double areaPBC = Vector.CrossProduct(PB, PC);
            double areaPCA = Vector.CrossProduct(PC, PA);

            double s = areaPBC / areaABC;
            double t = areaPCA / areaABC;
            double u = 1 - s - t;

            return new Tuple<double, double, double>(s, t, u);
        }
    }
}
