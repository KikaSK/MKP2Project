using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms;
using System.Windows.Forms.Integration;
using System.Windows.Media.Media3D;
using System.Diagnostics;

// pouzitie potrebnych kniznic OpenTK
using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

namespace MKP2___Template
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        /////////////////////////////////////////////////////////
        //                                                     //
        //                   GLOBALNE PREMENNE                 //
        //                                                     //
        /////////////////////////////////////////////////////////

        GLControl glControl;
       
        // Our two patches are store in the list "Patches"
        List<Patch> Patches = new List<Patch>();

        // The types of the left and of the right patch are stored in the variables "left" and "right". This will be useful if you want to change the types of the patch in user interface
        type left = type.TENSOR, right = type.TENSOR;

        // camera settings
        double Dist = new double(), Phi = new double(), Theta = new double(), oPhi = new double(), oTheta = new double(), prevPhi = new double(), prevTheta = new double(), prevDist = new double();

        // number of samples of the patch
        int nSamples = 20;
        // nubmer of control vertices in the direction m and n, respectively
        int nDeg;  // ----  the triangle patch requires the third direction

        // mouse settings
        double RightX, RightY;
        bool IsLeftDown, IsRightDown;
        int ActivePoint, ActivePatch;

        // keyboard settings
        bool IsZ = true, IsY = false, IsX = false;

        bool C1continuity = false;

        bool DrawCoordinateAxes;

        //-----------------------------------------------------------------------------------------------------------------------

        public MainWindow()
        {
            InitializeComponent();

            IsLeftDown = false;
            IsRightDown = false;

            // set value to true if you want to draw coordinate axes, false if not
            // x - red, y - green, z - blue
            DrawCoordinateAxes = false;
            T1checkbox.IsChecked = true;

            // initialize the parameters
            InitializeParams(left, right);
        }

        // initialization of paramters, when the application is launched or the patch is erased 
        private void InitializeParams(type Left, type Right)
        {
            Patches.Clear();
           
            //nSamples = 20;
            nDeg = Convert.ToInt32(Mbox.Text);

            // Defining the first patch 
            float[] color = { 0.2f, 0.6f, 0.7f };
            Patches.Add(new Patch(Left, nDeg, nSamples, color, placement.LEFT));

            // Defining the second patch
            float[] color2 = { 0.8f, 0.3f, 0.3f };
            Patches.Add(new Patch(Right, nDeg, nSamples, color2, placement.RIGHT));
            SetUp();
        }

        private void SetUp()
        {
            if (left == type.TRIANGLE && right == type.TENSOR)
            {
                for (int i = 0; i < nDeg + 1; ++i)
                {
                    Patches[0].ControlNet.Coordinates[i] =
                    Patches[1].ControlNet.Coordinates[i];
                }
                if (C1continuity)
                {
                    UpdateC1TensorTriangle();
                }
            }
        }

        // update TSBZ when TBZ changes
        private void UpdateC1TensorTriangle(int? index = null)
        {
            // subdivision - increaseng the degree
            Vector3 sub_point = new Vector3();
            sub_point = Patches[0].ControlNet.Coordinates[nDeg + 1];
            Patches[1].ControlNet.Coordinates[nDeg + 1] = 2 * Patches[0].ControlNet.Coordinates[0] - sub_point;
            for (int i = 1; i < nDeg + 1; ++i)
            {   
                    float coefA = (float)i / (float)nDeg;
                    float coefB = (float)(nDeg - i) / (float)nDeg;
                    Vector3 A = Patches[0].ControlNet.Coordinates[nDeg + i];
                if (coefB != 0)
                {
                    Vector3 B = Patches[0].ControlNet.Coordinates[nDeg + i + 1];
                    sub_point = coefA * A + coefB * B;
                }
                else
                    sub_point = A;
                    
                    // update TSBZ
                    Patches[1].ControlNet.Coordinates[i + nDeg + 1] = 2 * Patches[0].ControlNet.Coordinates[i] - sub_point;
                
            }
        }

        // Checking the compatibility
        public void CheckCompatibility()
        {
            if (IsLeftDown)
            {
                if (left == type.TENSOR && right == type.TENSOR)
                {
                    for (int i = 0; i < nDeg + 1; i++)
                    {
                        float heightDiff = Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[i]].Z
                            - Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[i]].Z;


                        if (Math.Abs(heightDiff) > float.Epsilon / 100)
                        {
                            Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[i]] =
                            Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[i]];

                            // move by the same vector
                            if (C1continuity)
                            { 
                                Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]] = new Vector3(
                                    Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].X,
                                    Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].Y,
                                    Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].Z - heightDiff);

                                Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]] = new Vector3(
                                    Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].X,
                                    Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].Y,
                                    Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].Z - heightDiff);
                            }


                        }
                        else if (C1continuity)
                        {
                            // changing patch 0
                            if (ActivePatch == 0)
                            {
                                heightDiff = Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].Z
                                   - Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[i]].Z;
                                // updating patch 1
                                Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]] = new Vector3(
                                    Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].X,
                                    Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].Y,
                                    Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[i]].Z - heightDiff
                                    );
                            }
                            // changing patch 1
                            else if (ActivePatch == 1)
                            {
                                heightDiff = Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]].Z
                                   - Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[i]].Z;
                                // updating patch 0
                                Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]] = new Vector3(
                                    Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].X,
                                    Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]].Y,
                                    Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[i]].Z - heightDiff
                                    );
                            }
                        }
                    }
                }
                else if (left == type.TRIANGLE && right == type.TRIANGLE)
                {
                    // moving common point 
                    if (ActivePoint <= nDeg)
                    {
                        float heightDiff = Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[ActivePoint]].Z
                             - Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[ActivePoint]].Z;

                        if (Math.Abs(heightDiff) > float.Epsilon / 100)
                        {
                            // C0 continuity
                            Patches[0].ControlNet.Coordinates[Patches[0].CommonEdge[ActivePoint]] =
                                Patches[1].ControlNet.Coordinates[Patches[1].CommonEdge[ActivePoint]];

                            // C1 continuity
                            if (C1continuity)
                            {
                                for (int i = ActivePoint; i < nDeg; ++i)
                                {
                                    Vector3 mid = (Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]] +
                                        Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]]) / 2;
                                    Patches[0].ControlNet.Coordinates[i + 1] = mid + (mid - Patches[0].ControlNet.Coordinates[i]);
                                    Patches[1].ControlNet.Coordinates[i + 1] = Patches[0].ControlNet.Coordinates[i + 1];
                                }
                                for (int i = ActivePoint - 1; i >= 0; i--)
                                {
                                    Vector3 mid = (Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[i]] +
                                        Patches[1].ControlNet.Coordinates[Patches[1].ParallelEdge[i]]) / 2;
                                    Patches[0].ControlNet.Coordinates[i] = mid + (mid - Patches[0].ControlNet.Coordinates[i + 1]);
                                    Patches[1].ControlNet.Coordinates[i] = Patches[0].ControlNet.Coordinates[i];
                                }
                            }
                        }
                    }
                    // moving parallel point
                    else if (ActivePoint <= 2 * nDeg && C1continuity)
                    {
                        int NonActivePatch = ActivePatch == 0 ? 1 : 0;
                        int index = ActivePoint - nDeg - 1;
                        Vector3 mid = (Patches[0].ControlNet.Coordinates[index] + Patches[0].ControlNet.Coordinates[index + 1]) / 2;
                        Patches[NonActivePatch].ControlNet.Coordinates[ActivePoint] = 2 * mid - Patches[ActivePatch].ControlNet.Coordinates[ActivePoint];
                    }
                }
                else if (left == type.TRIANGLE && right == type.TENSOR)
                {
                    // moving common point
                    if(ActivePoint < nDeg+1)
                    {
                        Patches[0].ControlNet.Coordinates[ActivePoint] =
                        Patches[1].ControlNet.Coordinates[ActivePoint];
                        if(C1continuity)
                        {
                            // subdivided point
                            Vector3 sub_point = new Vector3();
                            if (ActivePoint == 0) sub_point = Patches[0].ControlNet.Coordinates[nDeg+1];
                            else if (ActivePoint == nDeg) sub_point = Patches[0].ControlNet.Coordinates[2*nDeg];
                            else
                            {
                                float coefA = (float)ActivePoint / (float)nDeg;
                                float coefB = (float)(nDeg-ActivePoint) / (float)nDeg;
                                Vector3 A = Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[ActivePoint - 1]];
                                Vector3 B = Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[ActivePoint]];

                                sub_point = coefA * A + coefB * B;
                            }
                            // update TSBZ
                            Patches[1].ControlNet.Coordinates[ActivePoint+nDeg+1] = 2*Patches[0].ControlNet.Coordinates[ActivePoint]-sub_point;
                            
                        }
                    }
                    // moving parallel point on TBZ patch
                    else if(ActivePatch == 0 && ActivePoint <= 2*nDeg && C1continuity)
                    {
                        UpdateC1TensorTriangle();
                    }
                    // moving parallel point on TSBZ patch
                    else if (ActivePatch == 1 && ActivePoint <= 2*nDeg + 1 && C1continuity)
                    {
                        int index = ActivePoint - nDeg - 1;
                        // subdivided point s'
                        Vector3 sub_point = new Vector3();
                        sub_point = 2 * Patches[0].ControlNet.Coordinates[index] - Patches[1].ControlNet.Coordinates[ActivePoint];
                        // r_new
                        Vector3 r_new = new Vector3();
                        if(index == 0 || index == nDeg)
                        {
                            r_new = sub_point;
                        }
                        else
                        {
                            r_new = (nDeg * sub_point - index * Patches[0].ControlNet.Coordinates[nDeg + index]) / (nDeg - index);
                        }
                        if(index < nDeg)
                            Patches[0].ControlNet.Coordinates[nDeg + index + 1] = r_new;
                        else
                            Patches[0].ControlNet.Coordinates[nDeg + index] = r_new;
                       
                        UpdateC1TensorTriangle(index);
                    }

                }

                // --------------- !!! TODO !!! -------------------         
                //
                // Here insert your code for checking the continuity and adjusting the control
                // vertices to keep the desired degree of coninuity
                //
                // Example for vertex access:
                //      Say, we want to get the coordinates of the second control vertex of the
                //      second column of the patch which is on the left side and insert into the vector Output
                //      The code will look like following:
                // Vector3 Output = Patches[0].ControlNet.Coordinates[Patches[0].ParallelEdge[1]];   
                //
                // ------------------------------------------------
            }
        }

        private void Mplus_Click(object sender, RoutedEventArgs e)
        {
            nDeg++;
            Mbox.Text = Convert.ToString(nDeg);

            InitializeParams(left, right);
            // redraw the scene
            glControl.Invalidate();
        }


        private void Mminus_Click(object sender, RoutedEventArgs e)
        {
            if (nDeg > 0) nDeg--;
            Mbox.Text = Convert.ToString(nDeg);

            InitializeParams(left, right);
            // redraw the scene
            glControl.Invalidate();
        }


        private void Uminus_Click(object sender, RoutedEventArgs e)
        {
            if (nSamples > 0) nSamples--;
            Ubox.Text = Convert.ToString(nSamples);
            InitializeParams(left, right);
            glControl.Invalidate();
        }

        private void C0Checked(object sender, RoutedEventArgs e)
        {
            C1continuity = false;
            InitializeParams(left, right);
            glControl.Invalidate();
        }

        private void C1Checked(object sender, RoutedEventArgs e)
        {
            C1continuity = true;
            InitializeParams(left, right);
            glControl.Invalidate();
        }

        private void T1_Checked(object sender, RoutedEventArgs e)
        {
            left = type.TENSOR;
            right = type.TENSOR;
            T2checkbox.IsChecked = false;
            T3checkbox.IsChecked = false;
            InitializeParams(left, right);
            glControl.Invalidate();
        }
        private void T2_Checked(object sender, RoutedEventArgs e)
        {
            left = type.TRIANGLE;
            right = type.TENSOR;
            T1checkbox.IsChecked = false;
            T3checkbox.IsChecked = false;
            InitializeParams(left, right);
            glControl.Invalidate();
        }
        private void T3_Checked(object sender, RoutedEventArgs e)
        {
            left = type.TRIANGLE;
            right = type.TRIANGLE; 
            T2checkbox.IsChecked = false;
            T1checkbox.IsChecked = false;
            InitializeParams(left, right);
            glControl.Invalidate();
        }


        private void Uplus_Click(object sender, RoutedEventArgs e)
        {
            nSamples++;
            Ubox.Text = Convert.ToString(nSamples);
            InitializeParams(left, right);
            glControl.Invalidate();
        }
        

       






        //-----------------------------------------------------------------------------------------------------------------------

        /////////////////////////////////////////////////////////
        //                                                     //
        //                      PROCEDURY                      //
        //                                                     //
        /////////////////////////////////////////////////////////

        



        //-----------------------------------------------------------------------------------------------------------------------

        // draw the coordinate axes
        private void DrawAxes()
        {
            GL.Begin(PrimitiveType.Lines);
            GL.Color3(1.0f, 0.0f, 0.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(1.0f, 0.0f, 0.0f);
            GL.Vertex3(2.0f, 0.0f, 0.0f);

            GL.Color3(0.0f, 1.0f, 0.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(0.0f, 1.0f, 0.0f);
            GL.Vertex3(0.0f, 2.0f, 0.0f);

            GL.Color3(0.0f, 0.0f, 1.0f);
            GL.Vertex3(0.0f, 0.0f, 0.0f);
            GL.Color3(0.0f, 0.0f, 1.0f);
            GL.Vertex3(0.0f, 0.0f, 2.0f);
            GL.End();
        }


        //-----------------------------------------------------------------------------------------------------------------------
        //                                                      DRAWING
        //-----------------------------------------------------------------------------------------------------------------------
        
        // drawing the patch
        private void DrawPatch(Patch _patch)
        {
            // drawing triangles / quadrilaterals 
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill); // enabble filling of shapes with color 

            // color -- !!! TODO !!! -- edit if you want something different :-) 
            float[] diffuse = { 0.9f, 0.9f, 0.9f, 1.0f };
            float[] specular = { 0.1f, 0.1f, 0.1f, 0.5f };
            

            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Diffuse, diffuse);
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Specular, specular);
            
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Shininess, 0.1f);

            PrimitiveType prim = new PrimitiveType();
            if (_patch.TypeOfPatch == type.TENSOR) prim = PrimitiveType.Quads;
            if (_patch.TypeOfPatch == type.TRIANGLE) prim = PrimitiveType.Triangles;

            GL.Begin(prim); 
            for (int i = 0; i < _patch.Sampling.Indices.Count; i++)
            {
                
                GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Ambient, _patch.Color);
                
                GL.Normal3(0.0f, 0.0f, 1.0f);
                GL.Vertex3(_patch.Sampling.Coordinates[_patch.Sampling.Indices[i]]);
            }
            GL.End();

            // drawing the wireframe model
            GL.Translate(0.0f, 0.0f, 0.01f);
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line);

            float[] black = { 0.0f, 0.0f, 0.0f, 1.0f };
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Specular, black);
            GL.Material(MaterialFace.FrontAndBack, MaterialParameter.Ambient, black);
            
            GL.LineWidth(0.5f);

            GL.Begin(prim); 
            for (int i = 0; i < _patch.Sampling.Indices.Count; i++)
                GL.Vertex3(_patch.Sampling.Coordinates[_patch.Sampling.Indices[i]]);
            GL.End();
        }

//-----------------------------------------------------------------------------------------------------------------------

        // drawing the ControlNet
        private void DrawNet(Patch _patch)
        {
            // firstly, draw the wireframe of the control net
            GL.Translate(0.0f, 0.0f, 0.01f);
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Line); // zabezpeci vykreslenie drotoveho modelu
            
            GL.LineWidth(2.0f);
            GL.Color3(0.529f, 0.904f, 0.971f); // color of the wireframe net

            PrimitiveType prim = new PrimitiveType();
            if (_patch.TypeOfPatch == type.TENSOR) prim = PrimitiveType.Quads;
            if (_patch.TypeOfPatch == type.TRIANGLE) prim = PrimitiveType.Triangles;

            GL.Begin(prim); 
            for (int i = 0; i < _patch.ControlNet.Indices.Count; i++)
                GL.Vertex3(_patch.ControlNet.Coordinates[_patch.ControlNet.Indices[i]]);
            GL.End();

            // now, draw the control points of the control net
            GL.PointSize(6.0f);
            GL.Color3(0.490f, 0.116f, 0.116f); // color of the control points
            GL.Disable(EnableCap.DepthTest); // depth test is disabled, so the points are not  covered by the patch
            GL.Begin(PrimitiveType.Points);
            for (int i = 0; i < _patch.ControlNet.Coordinates.Count; i++)
                GL.Vertex3(_patch.ControlNet.Coordinates[i]);
            GL.End();
            GL.Enable(EnableCap.DepthTest); // enable depth test
        }

//-----------------------------------------------------------------------------------------------------------------------

        
//-----------------------------------------------------------------------------------------------------------------------

       

//-----------------------------------------------------------------------------------------------------------------------

        // drawing 
        private void GLControl_Paint(object sender, PaintEventArgs e)
        {
            // Modelview matrix
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            Matrix4 matLook = Matrix4.LookAt((float)(Dist * Math.Cos(Theta) * Math.Cos(Phi)), (float)(Dist * Math.Sin(Phi) * Math.Cos(Theta)), (float)(Dist * Math.Sin(Theta)), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
            GL.LoadMatrix(ref matLook);

            // perspective projection
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            Matrix4 matPers = Matrix4.CreatePerspectiveFieldOfView(0.785f, (float)glControl.Width / (float)glControl.Height, 0.1f, 10.5f);
            GL.LoadMatrix(ref matPers);

            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
            
            CheckCompatibility();
            if (DrawCoordinateAxes) DrawAxes();


            for (int i = 0; i < Patches.Count; i++)
                Patches[i].RecomputePatch();

            GL.Enable(EnableCap.Lighting);
            GL.Enable(EnableCap.DepthTest);
            for(int i = 0; i < Patches.Count; i++)
            DrawPatch(Patches[i]);
            GL.Disable(EnableCap.Lighting);
            for (int i = 0; i < Patches.Count; i++)
                DrawNet(Patches[i]);
            

            // the buffers need to swapped, so the scene is drawn
            glControl.SwapBuffers();
        }

//-----------------------------------------------------------------------------------------------------------------------

        // initialization of the window, where OpenTK drawing is used 
        private void WindowsFormsHost_Initialized(object sender, EventArgs e)
        {
            // Inicializacia OpenTK;
            OpenTK.Toolkit.Init();
            var flags = GraphicsContextFlags.Default;
            glControl = new GLControl(new GraphicsMode(32, 24), 2, 0, flags);
            glControl.MakeCurrent();
            glControl.Paint += GLControl_Paint;
            glControl.Dock = DockStyle.Fill;
            (sender as WindowsFormsHost).Child = glControl;

            // user controls
            glControl.MouseDown += GLControl_MouseDown;
            glControl.MouseMove += GLControl_MouseMove;
            glControl.MouseUp += GLControl_MouseUp;
            glControl.MouseWheel += GLControl_MouseWheel;

            // shading
            GL.ShadeModel(ShadingModel.Smooth);

            // color of the window
            GL.ClearColor(1.0f, 1.0f, 1.0f, 1.0f);

            
            GL.ClearDepth(1.0f);

            //enable z-buffering
            GL.Enable(EnableCap.DepthTest);
            GL.DepthFunc(DepthFunction.Lequal);
            GL.Hint( HintTarget.PerspectiveCorrectionHint, HintMode.Nicest);

            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            //smoothing
            GL.Enable(EnableCap.LineSmooth);
            GL.Enable(EnableCap.PointSmooth);

            // illumination
            float[] light_ambient = { 0.3f, 0.3f, 0.3f, 1.0f };
            float[] light_diffuse = { 0.4f, 0.4f, 0.4f, 0.0f };
            float[] light_specular = { 0.5f, 0.5f, 0.5f, 1.0f };
            float[] light_position = { 10.0f, 10.0f, 200.0f };
            GL.Light(LightName.Light0, LightParameter.Ambient, light_ambient);
            GL.Light(LightName.Light0, LightParameter.Diffuse, light_diffuse);
            GL.Light(LightName.Light0, LightParameter.Specular, light_specular);
            GL.Light(LightName.Light0, LightParameter.ConstantAttenuation, 1.0f);
            GL.Light(LightName.Light0, LightParameter.Position, light_position);
            GL.Enable(EnableCap.Light0);

            // parameters for the camera
            Phi = 0.6f; Theta = 0.6f; Dist = 3.8f;


        }

//-----------------------------------------------------------------------------------------------------------------------

        /////////////////////////////////////////////////////////
        //                                                     //
        //                 USER INTERFACE CONTROLS             //
        //                                                     //
        /////////////////////////////////////////////////////////
        private void GLControl_MouseDown(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Right) // camera is adjusted using RMB
            {
                IsRightDown = true;
                RightX = e.X;
                RightY = e.Y;
                oPhi = Phi;
                oTheta = Theta;
            }
            else if (e.Button == MouseButtons.Left) // using LMB we search for the control point beneath the mouse cursor 
            {
                // the idea of the searching -- when I am doing the inverse projection,
                // what points lie in the ray which is casted from the point beneath the cursor.
                // If there are any, I choose the closest one. 
                
                Vector3 start, end;

                int[] viewport = new int[4];
                Matrix4 modelMatrix, projMatrix;

                GL.GetFloat(GetPName.ModelviewMatrix, out modelMatrix);
                GL.GetFloat(GetPName.ProjectionMatrix, out projMatrix);
                GL.GetInteger(GetPName.Viewport, viewport);

                start = UnProject(new Vector3(e.X, e.Y, 0.0f), projMatrix, modelMatrix, new Size(viewport[2], viewport[3]));
                end = UnProject(new Vector3(e.X, e.Y, 1.0f), projMatrix, modelMatrix, new Size(viewport[2], viewport[3]));

                double se = Math.Sqrt(Vector3.Dot(start - end, start - end));
                for(int k = 0; k < Patches.Count; k++)
                for(int i = 0; i < Patches[k].ControlNet.Coordinates.Count; i++)
                {
                    double sA = Math.Sqrt(Vector3.Dot(Patches[k].ControlNet.Coordinates[i] - start, Patches[k].ControlNet.Coordinates[i] - start));
                    double eA = Math.Sqrt(Vector3.Dot(Patches[k].ControlNet.Coordinates[i] - end, Patches[k].ControlNet.Coordinates[i] - end));

                    if(sA + eA > se - 0.001 && sA + eA < se + 0.001)
                    {
                        ActivePoint = i;
                            ActivePatch = k;
                        IsLeftDown = true;

                        RightX = e.X;
                        RightY = e.Y;
                    }
                }
            }

            // redraw the scene
            glControl.Invalidate();
        }
        
        // Inverse projection
        public Vector3 UnProject(Vector3 mouse, Matrix4 projection, Matrix4 view, Size viewport)
        {
            Vector4 vec;

            vec.X = 2.0f * mouse.X / (float)viewport.Width - 1;
            vec.Y = -(2.0f * mouse.Y / (float)viewport.Height - 1);
            vec.Z = mouse.Z;
            vec.W = 1.0f;

            Matrix4 viewInv = Matrix4.Invert(view);
            Matrix4 projInv = Matrix4.Invert(projection);

            Vector4.Transform(ref vec, ref projInv, out vec);
            Vector4.Transform(ref vec, ref viewInv, out vec);

            if (vec.W > 0.000001f || vec.W < -0.000001f)
            {
                vec.X /= vec.W;
                vec.Y /= vec.W;
                vec.Z /= vec.W;
            }

            return vec.Xyz;
        }

//-----------------------------------------------------------------------------------------------------------------------

        private void GLControl_MouseMove(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (IsRightDown) // RMB - rotate the camera
            {
                IsRightDown = true;

                Phi = oPhi + (RightX - e.X) / 200.0f;
                Theta = oTheta + (e.Y - RightY) / 200.0f;
            }
            else if (IsLeftDown) // LMB - move the control vertex
            {
                IsLeftDown = true;

                float Scaling = 0.003f; 

                if (IsX)
                    Patches[ActivePatch].ControlNet.Coordinates[ActivePoint] = new Vector3(Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].X + Convert.ToSingle(RightX - e.X) * Scaling, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Y - Convert.ToSingle(RightY - e.Y) * Scaling, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Z);
                if (IsY)
                    Patches[ActivePatch].ControlNet.Coordinates[ActivePoint] = new Vector3(Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].X - Convert.ToSingle(RightY - e.Y) * Scaling, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Y - Convert.ToSingle(RightX - e.X) * Scaling, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Z);
                if (IsZ)
                    Patches[ActivePatch].ControlNet.Coordinates[ActivePoint] = new Vector3(Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].X, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Y, Patches[ActivePatch].ControlNet.Coordinates[ActivePoint].Z + Convert.ToSingle(RightY - e.Y) * Scaling);

                RightY = e.Y;
                RightX = e.X;
            }

            // redraw the scene
            glControl.Invalidate();
        }

//-----------------------------------------------------------------------------------------------------------------------

        private void GLControl_MouseUp(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Right) IsRightDown = false;
            if (e.Button == MouseButtons.Left) IsLeftDown = false;
        }

//-----------------------------------------------------------------------------------------------------------------------

        private void GLControl_MouseWheel(object sender, System.Windows.Forms.MouseEventArgs e)
        {
            Dist -= (double)e.Delta * 0.001; // zooming

            // redraw the scene
            glControl.Invalidate();
        }

//-----------------------------------------------------------------------------------------------------------------------

        private void Grid_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            GL.Viewport(0, 0, glControl.Width, glControl.Height);         
            
        }

        private void Window_KeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if (e.Key == Key.X) // view from above 1
            {
                if (IsX)
                {
                    IsX = false;
                    Phi = prevPhi;
                    Theta = prevTheta;
                    Dist = prevDist;

                    XLabelY.Visibility = System.Windows.Visibility.Hidden;
                    XLabelX.Visibility = System.Windows.Visibility.Hidden;
                    XRectY.Visibility = System.Windows.Visibility.Hidden;
                    XRectX.Visibility = System.Windows.Visibility.Hidden;
                }
                else
                {
                    IsX = true;
                    IsY = false;
                    IsZ = false;
                    prevPhi = Phi;
                    prevTheta = Theta;
                    prevDist = Dist;
                    Phi = 1.57;
                    Theta = 1.24;
                    Dist = 3.5;

                    LabelZ.Visibility = System.Windows.Visibility.Hidden;
                    RectZ.Visibility = System.Windows.Visibility.Hidden;
                    XLabelY.Visibility = System.Windows.Visibility.Visible;
                    XLabelX.Visibility = System.Windows.Visibility.Visible;
                    YLabelY.Visibility = System.Windows.Visibility.Hidden;
                    YLabelX.Visibility = System.Windows.Visibility.Hidden;
                    XRectY.Visibility = System.Windows.Visibility.Visible;
                    XRectX.Visibility = System.Windows.Visibility.Visible;
                    YRectY.Visibility = System.Windows.Visibility.Hidden;
                    YRectX.Visibility = System.Windows.Visibility.Hidden;
                }
            }

            if (e.Key == Key.Y) // view from above 2
            {
                if (IsY)
                {
                    IsY = false;
                    Phi = prevPhi;
                    Theta = prevTheta;
                    Dist = prevDist;

                    YLabelY.Visibility = System.Windows.Visibility.Hidden;
                    YLabelX.Visibility = System.Windows.Visibility.Hidden;
                    YRectY.Visibility = System.Windows.Visibility.Hidden;
                    YRectX.Visibility = System.Windows.Visibility.Hidden;
                }
                else
                {
                    IsY = true;
                    IsX = false;
                    IsZ = false;
                    prevPhi = Phi;
                    prevTheta = Theta;
                    prevDist = Dist;
                    Phi = 0;
                    Theta = 1.3;
                    Dist = 3.5;

                    LabelZ.Visibility = System.Windows.Visibility.Hidden;
                    RectZ.Visibility = System.Windows.Visibility.Hidden;
                    XLabelY.Visibility = System.Windows.Visibility.Hidden;
                    XLabelX.Visibility = System.Windows.Visibility.Hidden;
                    YLabelY.Visibility = System.Windows.Visibility.Visible;
                    YLabelX.Visibility = System.Windows.Visibility.Visible;
                    XRectY.Visibility = System.Windows.Visibility.Hidden;
                    XRectX.Visibility = System.Windows.Visibility.Hidden;
                    YRectY.Visibility = System.Windows.Visibility.Visible;
                    YRectX.Visibility = System.Windows.Visibility.Visible;
                }

            }

            if (e.Key == Key.Z)
            {
                if (IsZ)
                {
                    IsZ = false;

                    LabelZ.Visibility = System.Windows.Visibility.Hidden;
                    RectZ.Visibility = System.Windows.Visibility.Hidden;
                }
                else
                {
                    IsZ = true;
                    IsY = false;
                    IsX = false;
                    LabelZ.Visibility = System.Windows.Visibility.Visible;
                    RectZ.Visibility = System.Windows.Visibility.Visible;
                    XLabelY.Visibility = System.Windows.Visibility.Hidden;
                    XLabelX.Visibility = System.Windows.Visibility.Hidden;
                    YLabelY.Visibility = System.Windows.Visibility.Hidden;
                    YLabelX.Visibility = System.Windows.Visibility.Hidden;
                    XRectY.Visibility = System.Windows.Visibility.Hidden;
                    XRectX.Visibility = System.Windows.Visibility.Hidden;
                    YRectY.Visibility = System.Windows.Visibility.Hidden;
                    YRectX.Visibility = System.Windows.Visibility.Hidden;
                }
            }

            // redraw the scene
            glControl.Invalidate();
        }


    }


}
