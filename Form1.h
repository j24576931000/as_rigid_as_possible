//==========================================================================
//  This is the implementation of As-Rigid-As-Possible Shape Manipulation
//==========================================================================
#pragma once
#define BOOST_ALL_DYN_LINK
//#pragma comment(linker, "/SUBSYSTEM:CONSOLE /ENTRY:WinMainCRTStartup")
//#define _USE_MATH_DEFINES
//#define NOMINMAX
#include <windows.h>
#include <vector>
#include <cmath>
#include <time.h>

#include <iostream>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/glut.h>
#include "ArapInteractor_2.h"
#include "vavImage.h"
#include "TriangulationCgal.h"
//#include <stb_image.h>
#pragma region Application variables

vavImage	  *ImageEdge=NULL;	   //find contour
TriangulationCgal *Triangulate=NULL;	   //Delaunay triangulation	

TriMesh2D	  *test_1=NULL;		   	
ShapeView	  *ShapeView_Object=NULL;
ArapInteractor    *Arap=NULL;

int		  flag=-1;
int		  mouseX,mouseY;
int       tim=0;
time_t tt;
int next_t= 0 ;
int start=0;
int end_t = 0;
int play = 0;
struct Mode_Display {
	bool openImg;
	bool triangulation;
};

Mode_Display Current_Display;
#pragma endregion

namespace As_rigid_as_test {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			test_1 = new TriMesh2D;
			ImageEdge= new vavImage;
			

			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}

	protected: 

	private: HKOGLPanel::HKOGLPanelControl^  hkoglPanelControl1;

	private: System::Windows::Forms::Button^  button2;
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::Button^  record_button;
	private: System::Windows::Forms::Button^  play_button;
	private: System::Windows::Forms::Button^  stop_button;






	private: System::Windows::Forms::Timer^  timer1;
	private: System::Windows::Forms::Button^  reset;

	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			HKOGLPanel::HKCOGLPanelCameraSetting^  hkcoglPanelCameraSetting1 = (gcnew HKOGLPanel::HKCOGLPanelCameraSetting());
			HKOGLPanel::HKCOGLPanelPixelFormat^  hkcoglPanelPixelFormat1 = (gcnew HKOGLPanel::HKCOGLPanelPixelFormat());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			this->hkoglPanelControl1 = (gcnew HKOGLPanel::HKOGLPanelControl());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->record_button = (gcnew System::Windows::Forms::Button());
			this->play_button = (gcnew System::Windows::Forms::Button());
			this->stop_button = (gcnew System::Windows::Forms::Button());
			this->timer1 = (gcnew System::Windows::Forms::Timer(this->components));
			this->reset = (gcnew System::Windows::Forms::Button());
			this->SuspendLayout();
			// 
			// hkoglPanelControl1
			// 
			this->hkoglPanelControl1->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			hkcoglPanelCameraSetting1->Far = 1000;
			hkcoglPanelCameraSetting1->Fov = 45;
			hkcoglPanelCameraSetting1->Near = -1000;
			hkcoglPanelCameraSetting1->Type = HKOGLPanel::HKCOGLPanelCameraSetting::CAMERATYPE::ORTHOGRAPHIC;
			this->hkoglPanelControl1->Camera_Setting = hkcoglPanelCameraSetting1;
			this->hkoglPanelControl1->Location = System::Drawing::Point(136, 13);
			this->hkoglPanelControl1->Name = L"hkoglPanelControl1";
			hkcoglPanelPixelFormat1->Accumu_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Alpha_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Stencil_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			this->hkoglPanelControl1->Pixel_Format = hkcoglPanelPixelFormat1;
			this->hkoglPanelControl1->Size = System::Drawing::Size(734, 682);
			this->hkoglPanelControl1->TabIndex = 9;
			this->hkoglPanelControl1->Load += gcnew System::EventHandler(this, &Form1::hkoglPanelControl1_Load);
			this->hkoglPanelControl1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &Form1::hkoglPanelControl1_Paint);
			this->hkoglPanelControl1->KeyPress += gcnew System::Windows::Forms::KeyPressEventHandler(this, &Form1::hkoglPanelControl1_KeyPress);
			this->hkoglPanelControl1->MouseClick += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::hkoglPanelControl1_MouseClick);
			this->hkoglPanelControl1->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::hkoglPanelControl1_MouseDown);
			this->hkoglPanelControl1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::hkoglPanelControl1_MouseMove);
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(12, 13);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(118, 107);
			this->button2->TabIndex = 11;
			this->button2->Text = L"Open Img";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &Form1::button2_Click_1);
			// 
			// button3
			// 
			this->button3->Location = System::Drawing::Point(12, 135);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(118, 114);
			this->button3->TabIndex = 12;
			this->button3->Text = L"triangulation";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &Form1::Button3_Click);
			// 
			// record_button
			// 
			this->record_button->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"record_button.BackgroundImage")));
			this->record_button->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->record_button->Location = System::Drawing::Point(24, 279);
			this->record_button->Name = L"record_button";
			this->record_button->Size = System::Drawing::Size(70, 70);
			this->record_button->TabIndex = 13;
			this->record_button->UseVisualStyleBackColor = true;
			this->record_button->Click += gcnew System::EventHandler(this, &Form1::record_button_Click);
			// 
			// play_button
			// 
			this->play_button->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"play_button.BackgroundImage")));
			this->play_button->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->play_button->Location = System::Drawing::Point(24, 355);
			this->play_button->Name = L"play_button";
			this->play_button->Size = System::Drawing::Size(70, 70);
			this->play_button->TabIndex = 14;
			this->play_button->UseVisualStyleBackColor = true;
			this->play_button->Click += gcnew System::EventHandler(this, &Form1::play_button_Click);
			// 
			// stop_button
			// 
			this->stop_button->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"stop_button.BackgroundImage")));
			this->stop_button->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->stop_button->Location = System::Drawing::Point(24, 438);
			this->stop_button->Name = L"stop_button";
			this->stop_button->Size = System::Drawing::Size(70, 70);
			this->stop_button->TabIndex = 15;
			this->stop_button->UseVisualStyleBackColor = true;
			this->stop_button->Click += gcnew System::EventHandler(this, &Form1::stop_button_Click);
			// 
			// timer1
			// 
			this->timer1->Tick += gcnew System::EventHandler(this, &Form1::timer1_Tick);
			// 
			// reset
			// 
			this->reset->Location = System::Drawing::Point(24, 546);
			this->reset->Name = L"reset";
			this->reset->Size = System::Drawing::Size(75, 23);
			this->reset->TabIndex = 16;
			this->reset->Text = L"button1";
			this->reset->UseVisualStyleBackColor = true;
			this->reset->Click += gcnew System::EventHandler(this, &Form1::reset_Click);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->ClientSize = System::Drawing::Size(887, 707);
			this->Controls->Add(this->reset);
			this->Controls->Add(this->stop_button);
			this->Controls->Add(this->play_button);
			this->Controls->Add(this->record_button);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->hkoglPanelControl1);
			this->Name = L"Form1";
			this->Text = L"As rigid as possible shape manipulation";
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->ResumeLayout(false);

		}
#pragma endregion
	private: System::Void Form1_Load(System::Object^  sender, System::EventArgs^  e) {
		 }
	private: System::Void dotNetBarManager1_ItemClick(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button1_MouseClick(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		 hkoglPanelControl1->Invalidate();
	 }
private: System::Void hkoglPanelControl1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {
		 Render_Init(hkoglPanelControl1->Size.Width, hkoglPanelControl1->Size.Height);
		 if(Current_Display.openImg)
		 {
			 ImageEdge->drawImage();
			 //ImageEdge->drawImage_rgb();
		 }
		 if(Current_Display.triangulation)
		 {
			 //ImageEdge->drawImage();
			 Arap->OnDraw(); 
			
		 }
		 int now = System::DateTime::Now.Second;
		 //std::cout << now << std::endl;
		 if (start == 1 && end_t == 0)
		 {
			 //std::cout << now << std::endl;
			 Arap->save_mesh();
			 //start = 0;
		 }
		 if (play == 1)
		 {
			 Arap->return_mesh(tim);
			 tim += 1;
		 }
		 if (play == 1 && tim == (Arap->record_mesh.size() - 1))
		 {
			 tim = 0;
		 }
		 //if (next_t == now)
		 //{
			 //start = 1;
		 //}
		  //next_t = now + 1;
		 hkoglPanelControl1->Invalidate();
	 }
private: System::Void hkoglPanelControl1_Load(System::Object^  sender, System::EventArgs^  e) {
		 OpenGLinitial();
		 PanelResize(hkoglPanelControl1->Size.Width , hkoglPanelControl1->Size.Height);
	 }
private: System::Void hkoglPanelControl1_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		 if(e->Button==System::Windows::Forms::MouseButtons::Left)//move the control point
		 {
			 if(flag!=-1)
			 {
				 //std::cout<<"position X:"<<e->X<<" Y:"<<e->Y<<std::endl;
				 //==========================
// 				 clock_t start, finish;
// 				 start = clock();
				 //==========================
				 Arap->OnMotion(e->X-50,e->Y-50,flag,1,1);//0.008s

				 //==============================
// 				 finish = clock();
// 				 cout <<"Time Consume :" <<(double)(finish - start) / CLOCKS_PER_SEC<<endl;
				 //==============================
				 hkoglPanelControl1->Invalidate();
			 }		
		 }

	 }
private: System::Void hkoglPanelControl1_MouseClick(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		 if(e->Button==System::Windows::Forms::MouseButtons::Left && mouseX==e->X && mouseY==e->Y)//add the control point
		 {	
			 Arap->OnMouse(0,1,e->X-50,e->Y-50);
		 }
		 if(e->Button==System::Windows::Forms::MouseButtons::Right)				  //delete the control point
		 {
			 Arap->OnMouse(2,1,e->X-50,e->Y-50);
		 }
		 hkoglPanelControl1->Invalidate();
	 }
private: System::Void hkoglPanelControl1_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		 if(e->Button==System::Windows::Forms::MouseButtons::Left)
		 {
			mouseX=e->X;
			mouseY=e->Y;
			flag=Arap->getVertex(e->X-50,e->Y-50);						 //get control point ID
		 }
	 }
private: System::Void button2_Click_1(System::Object^  sender, System::EventArgs^  e) {//open image 
		 Current_Display.openImg=1;
		 ImageEdge->ReadImage("eevee_mask.png");
		 //ImageEdge->ReadImage_rgb("1111.png");
		 *ImageEdge = (ImageEdge->CannyEdge());
	 
		 //std::cout<<"Load img :"<<ImageEdge->GetHeight()<<"*"<<ImageEdge->GetWidth()<<std::endl;
	 
		 hkoglPanelControl1->Invalidate();
	 }
private: System::Void Button3_Click(System::Object^  sender, System::EventArgs^  e) {//CGAL Delaunay Triangulation
		 if(Current_Display.openImg)
		 {
			 Current_Display.triangulation=1;
			 Current_Display.openImg=0;
			 Triangulate=new TriangulationCgal;
			 Triangles Tris;
			 Vector2s meshPointset;
			 Vector2s ContourPoint=ImageEdge->GetContour();

			 for(int i=0;i<ContourPoint.size();i+=15)
			 {
				  Triangulate->AddPoint(ContourPoint[i][0],ContourPoint[i][1]);
			 }

			 Triangulate->DelaunayMesher2();
			 meshPointset=Triangulate->MeshPointSet();
			 //std::cout << "meshPointset.size(): " << meshPointset.size() << std::endl;
			 for(int i=0 ;i < meshPointset.size() ; i++ )
			 {
				 test_1->vertices.push_back(Point2D(meshPointset[i][0],meshPointset[i][1]));
			 }
			 for (int i = 0; i < meshPointset.size(); i++)
			 {
				 test_1->texcoord.push_back(Point2D((meshPointset[i][0]-1)/535, (meshPointset[i][1]-1)/533));
			 }
 			 Tris=Triangulate->GetTriangles();
 			 //std::cout<<"Tris.size() :"<<Tris.size()<<std::endl;
			 Tri v;
			 for(int i=0;i<Tris.size();i++)
			 {
				 if(!ImageEdge->IsinsidePoint(Tris[Tris.size()-1-i].m_Points[0][0],Tris[Tris.size()-1-i].m_Points[0][1],
					 Tris[Tris.size()-1-i].m_Points[1][0],Tris[Tris.size()-1-i].m_Points[1][1],
					 Tris[Tris.size()-1-i].m_Points[2][0],Tris[Tris.size()-1-i].m_Points[2][1]))
 					 continue;
				 for(int j=0;j<3;j++)
				 {
					 v[j]=Triangulate->getVertexID(Tris[Tris.size()-1-i].m_Points[j][0],Tris[Tris.size()-1-i].m_Points[j][1]);
				 }
				test_1->tris.push_back(v);
			 }
 			 Arap= new ArapInteractor(ShapeView_Object,*test_1);

		 }
		 hkoglPanelControl1->Invalidate();
	 }
		private: System::Void record_button_Click(System::Object^  sender, System::EventArgs^  e) {
			//start = clock();
			//tt = time(NULL);
			//this->timer1->Enabled = true;

			//this->timer1->Enabled = false;
			if (start == 0&& end_t == 0)
			{
				//this->timer1->Enabled = true;
				tim = 0;
				start = 1;
				end_t = 0;
			}
			else if (start == 1 && end_t == 0)
			{
				start = 0;
				end_t = 1;
			}
				
			hkoglPanelControl1->Invalidate();
		}
		
		private: System::Void play_button_Click(System::Object^  sender, System::EventArgs^  e) {
			play = 1;
			//this->timer1->Enabled = true;
			hkoglPanelControl1->Invalidate();
		}
		private: System::Void stop_button_Click(System::Object^  sender, System::EventArgs^  e) {
			play =0;
			hkoglPanelControl1->Invalidate();
		}
		private: System::Void timer1_Tick(System::Object^  sender, System::EventArgs^  e) {

			//DateTime now = DateTime.Now;
			//int now = System::DateTime::Now.Second;
			////std::cout << now << std::endl;
			//if (start == 1&&end_t == 0)
			//{			
			//	//std::cout << now << std::endl;
			//	Arap->save_mesh();
			//	//start = 0;
			//}
			//if (play == 1)
			//{	
			//	Arap->return_mesh(tim);
			//	tim += 1;
			//}
			//if (play == 1&&tim==(Arap->record_mesh.size()-1) )
			//{
			//	tim = 0;
			//}
			////if (next_t == now)
			////{
			//	//start = 1;
			////}
			// //next_t = now + 1;
			//hkoglPanelControl1->Invalidate();
			
		}
		private: System::Void reset_Click(System::Object^  sender, System::EventArgs^  e) {
			start = 0;
			end_t = 0;
			play = 0;
			Arap->recover_mesh();
			Arap->record_mesh.clear();
			
			hkoglPanelControl1->Invalidate();
		}
	private: System::Void hkoglPanelControl1_KeyPress(System::Object^  sender, System::Windows::Forms::KeyPressEventArgs^  e) {
	
		if (e->KeyChar == 'f' || e->KeyChar == 'F') {
			
			Arap->OnKeyboard(e->KeyChar);
			hkoglPanelControl1->Invalidate();
			//e.Handled = true;
		}
		else if (e->KeyChar == '1' ) {
			Arap->OnKeyboard(e->KeyChar);

			hkoglPanelControl1->Invalidate();
			//e.Handled = true;
		}
	}
};
}

