using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK;



namespace MKP2___Template
{
    // enum type  determines the type of a patch
    public enum type
    {
        TENSOR,  // Bezier tensor product patch
        TRIANGLE // Bezier triangle patch
    };

    // enum placement determines the position of a patch
    public enum placement 
    {
        LEFT,
        MIDDLE,
        RIGHT
    }

    class Patch
    {
    
        // Class CNet describes the control vertices of a patch
        public class CNet
        {
            public List<int> Indices;
            public List<Vector3> Coordinates;

            public CNet()
            {
                Indices = new List<int>(); // Indices of a control net
                Coordinates = new List<Vector3>(); // Coordinates of control vertices
            }
        }

        // sampled patch (points computed using the de Casteljau algorithm)
        public class Sampl
        {
            public List<int> Indices; 
            public List<Vector3> Coordinates;

            public Sampl()
            {
                Indices = new List<int>();
                Coordinates = new List<Vector3>();
            }
        }

        public type TypeOfPatch; 
        public placement Place;
        public int NumberOfSamples, Degree;

        public CNet ControlNet;
        public Sampl Sampling;
        public float[] Color;
        public List<int> CommonEdge, ParallelEdge; // CommonEdge are the indices of those vertices which are "common" for both patches, ParallelEdge are the indices of those vertices which are "in the second column" of a patch

        // Initialization of a patch
        public Patch(type _TypeOfPatch, int _Degree, int _NumberOfSamples, float[] _Color, placement _Place)
        {
            TypeOfPatch = _TypeOfPatch;
            NumberOfSamples = _NumberOfSamples;
            Degree = _Degree;
            Color = _Color;
            Place = _Place;

            ControlNet = new CNet();
            Sampling = new Sampl();
            CommonEdge = new List<int>();
            ParallelEdge = new List<int>();

            // Initial sampling of a patch
            Sampling.Coordinates = Sample(TypeOfPatch, NumberOfSamples, NumberOfSamples);
            Sampling.Indices = GetIndices(TypeOfPatch, NumberOfSamples, NumberOfSamples, true);

            // Initial control net
            ControlNet.Coordinates = Sample(TypeOfPatch, Degree, Degree);
            ControlNet.Indices = GetIndices(TypeOfPatch, Degree, Degree, false);

            // --------------- !!! TODO !!! -------------------
            //
            // According to the value "Place" embed the patch to the correct place
            // in the scene. You can transform the whole control net using the
            // affine transformations.
            // Also, fill the lists "CommonEdge" and "ParallelEdge" with proper
            // indices of the corresponding vertices (take into account "TypeOfPatch")
            //
            // ------------------------------------------------

            float kScale = 1.0f;
            float kTranslate = -0.5f;
            if (Place == placement.LEFT)
            {
                if (TypeOfPatch == type.TENSOR)
                {
                    Scale(kScale, kScale, 0);
                    Translate(kScale + kTranslate, 0, 0);
                    //Console.WriteLine("Tensor");
                    for (int i = 0; i<Degree+1; i++)
                    {
                        //Console.WriteLine(i + " " + ControlNet.Coordinates[i]);
                        CommonEdge.Add(i);
                        ParallelEdge.Add(Degree + 1 + i);
                    }
                }
                else
                {
                    RotateZ(-(float)Math.PI / 2);
                    Scale(kScale, -kScale, 0);
                    Translate(-kTranslate, kScale / (Degree), 0);
                    
                    for (int i = 0; i<Degree; i++)
                    {
                        //Console.WriteLine(i + " " + ControlNet.Coordinates[i]);
                        CommonEdge.Add(i);
                        ParallelEdge.Add(i + Degree + 1);
                        //Console.WriteLine(i + " parallel " + ControlNet.Coordinates[i + Degree + 1]);
                    }
                    Console.WriteLine(Degree + " " + ControlNet.Coordinates[Degree]);
                    CommonEdge.Add(Degree);
                }
            }
            else if (Place == placement.RIGHT)
            {
                if (TypeOfPatch == type.TENSOR)
                {
                    Scale(-kScale, kScale, 0);
                    Translate(-kScale + kTranslate, 0, 0);

                    //Console.WriteLine("Tensor");
                    for (int i = 0; i < Degree+1; i++)
                    {
                        //Console.WriteLine(i + " " + ControlNet.Coordinates[Degree + 1 + i]);
                        CommonEdge.Add(i);
                        ParallelEdge.Add(Degree + 1 + i);
                    }
                }
                else
                {
                    RotateZ((float)Math.PI / 2);
                    Scale(kScale, kScale, 0);
                    Translate(3*kTranslate, kScale / (Degree), 0);
                    //TranlateTriangleTenzor(kScale / (Degree));

                    //Console.WriteLine("TriangleR");
                    for (int i = 0; i < Degree; i++)
                    {
                        //Console.WriteLine(i + " " + ControlNet.Coordinates[i]);
                        CommonEdge.Add(i);
                        ParallelEdge.Add(i + Degree + 1);
                        //Console.WriteLine(i + " parallel " + ControlNet.Coordinates[i + Degree + 1]);
                    }
                    //Console.WriteLine(Degree + " " + ControlNet.Coordinates[Degree]);
                    CommonEdge.Add(Degree);
                }
            }

        }

        // sampling of the initial patch with the given number of samples eU, eV
        // in the direction u, v, respectively
        private List<Vector3> Sample(type _TypeOfPatch, int eU, int eV)
        {
            List<Vector3> SampleList = new List<Vector3>();

            if (_TypeOfPatch == type.TENSOR)
            {
                for (int i = 0; i <= eV; i++)
                    for (int j = 0; j <= eU; j++)
                        SampleList.Add(new Vector3(-1.0f + 2.0f * i / eV, -1.0f + 2.0f * j / eU, 0.0f));
            }  
            
            if (_TypeOfPatch == type.TRIANGLE)
            {
                // --------------- !!! TODO !!! -------------------
                // 
                // Insert the sampling of the triangle from the previous assignment
                //
                // ------------------------------------------------
                Vector3 A = new Vector3(-1, -1, 0);
                Vector3 B = new Vector3(1, -1, 0);
                Vector3 C = new Vector3(0, 1, 0);

                for (int i = 0; i <= eU; ++i)
                {
                    // Points on the edges in the ith storey of the triangle pyramid
                    Vector3 P = A + 1.0f * i / eU * (C - A);
                    Vector3 Q = B + 1.0f * i / eU * (C - B);

                    for (int j = 0; j < (eU - i); ++j)
                    {
                        Vector3 R = P + 1.0f * j / (eU - i) * (Q - P);
                        SampleList.Add(R);
                    }
                    SampleList.Add(Q);
                }

                return SampleList;
            }

            return SampleList;
        }

        // getting indicies for a patch with given number of samples eU, eV in the direction u, v, respectively 
        private List<int> GetIndices(type _TypeOfPatch, int eU, int eV, bool DrawAll)
        {
            List<int> IndList = new List<int>();
            if (_TypeOfPatch == type.TENSOR)
            {
                // indices for rectangles - quadruples  
                for (int i = 0; i < eV; i++)
                    for (int j = 0; j < eU; j++)
                    {
                        IndList.Add(i * (eU + 1) + j);
                        IndList.Add(i * (eU + 1) + j + 1);
                        IndList.Add((i + 1) * (eU + 1) + j + 1);
                        IndList.Add((i + 1) * (eU + 1) + j);
                    }
            }

            if (_TypeOfPatch == type.TRIANGLE)
            {
                // --------------- !!! TODO !!! -------------------
                //
                //   insert the indices for the triangle patch from previous homework
                //
                // ------------------------------------------------
                int pos = 0;
                for (int i = 0; i < eU; ++i) // rows
                {
                    for (int j = 0; j < eU - i; ++j) // cols
                    {
                        IndList.Add(pos);
                        IndList.Add(pos + 1);
                        IndList.Add(pos + eU - i + 1);
                        pos++;
                    }
                    pos++;
                }

                if (DrawAll)
                {
                    pos = 0;
                    for (int i = 0; i < eU - 1; ++i)
                    {
                        for (int j = 0; j < eU - i - 1; ++j)
                        {
                            IndList.Add(pos + 1);
                            IndList.Add(pos + eU + 2 - i);
                            IndList.Add(pos + eU + 1 - i);
                            pos++;
                        }
                        pos += 2;
                    }
                }
            }

            return IndList;
        }

        // Computation of points of the patch 
        public void RecomputePatch()
        {
            // --------------- !!! TODO !!! -------------------
            //
            //   Uncomment the line below after filling the code for computation of the sampling
            //
            // ------------------------------------------------
            
            Sampling.Coordinates.Clear();
            Sampling.Coordinates = Sample(TypeOfPatch, NumberOfSamples, NumberOfSamples);

            if (TypeOfPatch == type.TENSOR)
            {
                // --------------- !!! TODO !!! -------------------
                //
                //   insert the computation of the tensor patch from previous homework
                //
                // ------------------------------------------------

                //Sampling.Normals.Clear();
                MyFunctions.GetSamplingAndNormals(ref Sampling.Coordinates, ref ControlNet.Coordinates, Degree + 1, Degree + 1);

            }


            if (TypeOfPatch == type.TRIANGLE)
            {
                // --------------- !!! TODO !!! -------------------
                //
                //   insert the computation of the triangle patch from previous homework
                //
                // ------------------------------------------------
                for (int i = 0; i < Sampling.Coordinates.Count; ++i)
                {
                    float s, t, u;
                    Tuple<double, double, double> BC = MyFunctions.GetBC(Sampling.Coordinates[i], new Vector3(-1, -1, 0), new Vector3(1, -1, 0), new Vector3(0, 1, 0));
                    s = (float)BC.Item1;
                    t = (float)BC.Item2;
                    u = (float)BC.Item3;

                    int n = Degree;
                    List<Vector3> CP = new List<Vector3>(ControlNet.Coordinates);
                    while (n >= 1)
                    {
                        List<Vector3> newCP = new List<Vector3>();
                        List<int> indices = GetIndices(TypeOfPatch, n, n, true);
                        for (int j = 0; j < indices.Count; j += 3)
                        {
                            Vector3 A = CP[indices[j]];
                            Vector3 B = CP[indices[j + 1]];
                            Vector3 C = CP[indices[j + 2]];
                            newCP.Add(s * A + t * B + u * C);
                        }
                        CP = new List<Vector3>(newCP);
                        newCP.Clear();
                        n--;
                        /*
                        if (n == 1)
                        {
                            Vector3 a = CP[1] - CP[0];
                            Vector3 b = CP[2] - CP[0];
                            Vector3 normal = (Vector3.Cross(a, b)).Normalized();
                            Sampling.Normals[i] = normal;
                        }
                        */
                    }
                    Sampling.Coordinates[i] = CP[0];
                }
            }
        }        
        
        public void Translate(float eX, float eY, float eZ)
        {
            for(int i = 0; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X + eX,
                                                          ControlNet.Coordinates[i].Y + eY,
                                                          ControlNet.Coordinates[i].Z + eZ);
            }
            for(int i = 0; i < Sampling.Coordinates.Count; i++)
            {
                Sampling.Coordinates[i] = new Vector3(Sampling.Coordinates[i].X + eX,
                                                        Sampling.Coordinates[i].Y + eY,
                                                        Sampling.Coordinates[i].Z + eZ);
            }
        }

        public void Scale(float eX, float eY, float eZ)
        {
            for (int i = 0; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X * eX,
                                                          ControlNet.Coordinates[i].Y * eY,
                                                          ControlNet.Coordinates[i].Z * eZ);
            }
            for (int i = 0; i < Sampling.Coordinates.Count; i++)
            {
                Sampling.Coordinates[i] = new Vector3(Sampling.Coordinates[i].X * eX,
                                                        Sampling.Coordinates[i].Y * eY,
                                                        Sampling.Coordinates[i].Z * eZ);
            }
        }

        public void RotateZ(float alpha)
        {
            //  cos alpha -sin alpha 0   
            //  sin alpha cost alpha 0
            //     0          0      1
            for (int i = 0; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(
                    ControlNet.Coordinates[i].X * (float)Math.Cos(alpha) - ControlNet.Coordinates[i].Y * (float)Math.Sin(alpha),
                    ControlNet.Coordinates[i].X * (float)Math.Sin(alpha) + ControlNet.Coordinates[i].Y * (float)Math.Cos(alpha),
                    ControlNet.Coordinates[i].Z);
            }
        }
        public void TranlateTriangleTenzor(float eY)
        {
            /*
            for(int i = 0; i<Degree; ++i)
            {
                ControlNet.Coordinates[i] = 
            }
            */
            for (int i = Degree+1; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X,
                                                          ControlNet.Coordinates[i].Y+eY,
                                                          ControlNet.Coordinates[i].Z);
            }
        }
        
    }
}
