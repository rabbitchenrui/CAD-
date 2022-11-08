#include <GLTools.h>
#include <stdio.h>
#include <gl/glew.h>
#include <gl/glut.h>
#include <math.h>
#include <iostream>
#include<vector>
using namespace std;
struct Solid;
struct Face;
struct Loop;
struct HalfEdge;
struct Vertex;
struct Edge;//声明使用到的元素
struct Solid
{
	int id;
	Face* faces; 
	Edge* edges; 
	Solid* next;
	Solid* pre;

	int vnum;//点的数量
	int fnum;//面的数量
	int lnum;//环的数量
	Solid() : id(0), faces(NULL), edges(NULL), next(NULL), pre(NULL), fnum(0), vnum(0), lnum(0) {}
};
struct Face
{
	int id;
	Solid* solid; 
	Loop* out_lp; 
	Loop* inner_lp;
	Face* next;
	Face* pre;
	int innum;//内环数量 扫掠用
	Face() : id(0), solid(NULL), out_lp(NULL), next(NULL), pre(NULL), inner_lp(NULL), innum(0) {}
};
struct Loop
{
	int id;
	HalfEdge* halfedges; 
	Face* face; 
	Loop* next;
	Loop* pre;

	Loop() : id(0), halfedges(NULL), face(NULL), next(NULL), pre(NULL) {}
};
struct Edge
{
	HalfEdge* half_l; 
	HalfEdge* half_r; 
	Edge* next;
	Edge* pre;
	Edge() : half_l(NULL), half_r(NULL), next(NULL), pre(NULL) {}
};
struct HalfEdge
{
	Edge* edge; 
	Vertex* sv; 
	Vertex* ev; 
	Loop* lp; 
	HalfEdge* next;
	HalfEdge* pre;
	HalfEdge* brother;
	HalfEdge() : edge(NULL), sv(NULL), lp(NULL), next(NULL), pre(NULL), brother(NULL) {}
};
struct Vertex
{
	int id;
	double coordinate[3];
	Vertex* next;
	Vertex* pre;
	Vertex(double x, double y, double z) : id(0), next(NULL), pre(NULL)
	{
		coordinate[0] = x;
		coordinate[1] = y;
		coordinate[2] = z;
	}
};//定义使用到的拓扑元素
void addEdgeIntoSolid(Edge* edge, Solid*& solid);
void addFaceIntoSolid(Face* face, Solid*& solid);
void addLoopIntoFace(Loop* loop, Face* face);//声明拓扑元素之间需要的操作 
vector<Vertex*> v_list;
vector<Loop*> l_list;
vector<Face*> sweep_list;//定义三个向量 分别保存 点/环/面 用于画图和扫掠
Solid* mvfs(double point[3], Vertex*& vertex);
HalfEdge* mev(Vertex* sv, double point[3], Loop* lp);
Loop* mef(Vertex* sv, Vertex* ev, Loop* lp, bool mark);
Loop* kemr(Vertex* sv, Vertex* ev, Loop* lp);
void kfmrh(Face* fa, Face* fb);//声明基本的欧拉操作
void sweep(double dir[3], double dist);//声明扫掠操作
Solid* mvfs(double point[3], Vertex*& vertex)
{
	Solid* solid = new Solid();
	Face* face = new Face();
	Loop* out_lp = new Loop();
	vertex = new Vertex(point[0], point[1], point[2]);

	vertex->id = solid->vnum;
	out_lp->id = solid->lnum;
	face->id = solid->fnum;

	l_list.push_back(out_lp);
	v_list.push_back(vertex);
	solid->vnum += 1;
	solid->fnum += 1;
	solid->lnum += 1;


	solid->faces = face;
	face->solid = solid;

	face->out_lp = out_lp;
	out_lp->face = face;

	return solid;
}
HalfEdge* mev(Vertex* sv, double point[3], Loop* loop)
{
	Solid* solid = loop->face->solid;
	Edge* edge = new Edge();
	HalfEdge* half_l = new HalfEdge();
	HalfEdge* half_r = new HalfEdge();
	Vertex* ev = new Vertex(point[0], point[1], point[2]);

	ev->id = solid->vnum;
	v_list.push_back(ev);
	solid->vnum += 1;

	half_l->sv = sv;
	half_l->ev = ev;
	half_r->sv = ev;
	half_r->ev = sv;

	edge->half_l = half_l;
	edge->half_r = half_r;
	half_l->edge = edge;
	half_r->edge = edge;

	half_r->brother = half_l;
	half_l->brother = half_r;

	half_l->lp = loop;
	half_r->lp = loop;


	if (loop->halfedges == NULL)
	{
		half_l->next = half_r;
		half_r->next = half_l;

		half_l->pre = half_r;
		half_r->pre = half_l;
		loop->halfedges = half_l;
	}
	else
	{
		HalfEdge* thalf = loop->halfedges;
		while (thalf->ev != sv)thalf = thalf->next;
		half_r->next = thalf->next;
		thalf->next->pre = half_r;
		thalf->next = half_l;
		half_l->pre = thalf;
		half_l->next = half_r;
		half_r->pre = half_l;
	}


	addEdgeIntoSolid(edge, solid);
	return half_l;
}
Loop* mef(Vertex* sv, Vertex* ev, Loop* loop, bool mark)
{
	Solid* solid = loop->face->solid;
	Edge* edge = new Edge();
	HalfEdge* half_l = new HalfEdge();
	HalfEdge* half_r = new HalfEdge();
	Loop* newLoop = new Loop();

	half_l->sv = sv;
	half_l->ev = ev;
	half_r->sv = ev;
	half_r->ev = sv;

	half_r->brother = half_l;
	half_l->brother = half_r;

	half_l->edge = edge;
	half_r->edge = edge;
	edge->half_l = half_l;
	edge->half_r = half_r;



	HalfEdge* thalf = loop->halfedges;
	HalfEdge* tmpa, * tmpb, * tmpc;
	while (thalf->ev != sv)thalf = thalf->next;
	tmpa = thalf;

	while (thalf->ev != ev)thalf = thalf->next;
	tmpb = thalf;

	thalf = thalf->next;
	while (thalf->ev != ev)thalf = thalf->next;
	tmpc = thalf;

//分环
	half_r->next = tmpa->next;
	tmpa->next->pre = half_r;
	tmpa->next = half_l;
	half_l->pre = tmpa;

	half_l->next = tmpb->next;
	tmpb->next->pre = half_l;
	tmpb->next = half_r;
	half_r->pre = tmpb;
	loop->halfedges = half_l;
	while (loop->halfedges->sv != ev)
	{
		loop->halfedges = loop->halfedges->next;
	}
	newLoop->halfedges = half_r;
	while (newLoop->halfedges->sv != ev)
	{
		newLoop->halfedges = newLoop->halfedges->next;
	}
	half_l->lp = loop;
	half_r->lp = newLoop;

	Face* face = new Face();

	newLoop->id = solid->lnum;
	solid->lnum += 1;
	l_list.push_back(newLoop);

	addFaceIntoSolid(face, solid);

	addLoopIntoFace(newLoop, face);

	if (tmpc == tmpb)
	{
		if (mark)//底面
		{
			sweep_list.push_back(half_l->lp->face);
		}
	}
	else
	{
		sweep_list.push_back(half_r->lp->face);
	}


	addEdgeIntoSolid(edge, solid);

	return loop;
}
Loop* kemr(Vertex* sv, Vertex* ev, Loop* loop)//sv must belong to the outer loop
{
	HalfEdge* tmpa, * tmpb, * hal;
	Face* face = loop->face;
	Loop* inlp = new Loop();
	Solid* solid = loop->face->solid;

	hal = loop->halfedges;

	while (hal->sv != sv || hal->ev != ev)hal = hal->next;
	tmpa = hal;

	while (hal->sv != ev || hal->ev != sv)hal = hal->next;
	tmpb = hal;

	tmpb->pre->next = tmpa->next;
	tmpa->pre->next = tmpb->next;

	loop->face->solid->faces->out_lp->halfedges = tmpa->pre;

	inlp->halfedges = tmpb->pre;
	tmpb->pre->lp = inlp;

	inlp->id = solid->lnum;
	solid->lnum += 1;
	l_list.push_back(inlp);

	addLoopIntoFace(inlp, tmpa->pre->brother->lp->face);

	delete tmpa;
	delete tmpb;

	return inlp;
}
void kfmrh(Face* fa, Face* fb)//fa indicate the outface, fb indicate the innerface
{
	Loop* loop = fb->out_lp;
	fa->solid->fnum -= 1;

	Solid* solid = fa->solid;
	Face* face = solid->faces;
	if (face == fb)
	{
		solid->faces = face->next;
	}
	else
	{
		Face* tf = face;
		while (face != fb && face != NULL)
		{
			tf = face;
			face = face->next;
		}
		tf->next = face->next;
	}
}
void sweep(double dir[3], double d)
{
	Vertex* startv, * nextv, * upv, * upprev;
	HalfEdge* he, * suphe, * uphe;
	double point[3];

	for (int i=0;i<sweep_list.size();i++)
	{

		//第一个点
		Loop* loop = (sweep_list[i])->out_lp;
		he = loop->halfedges;
		startv = he->sv;
		point[0] = startv->coordinate[0] + d * dir[0];
		point[1] = startv->coordinate[1] + d * dir[1];
		point[2] = startv->coordinate[2] + d * dir[2];

		suphe = mev(startv, point, loop);//构建第一个竖直边
		upprev = suphe->ev;//记录第一个顶面点
		he = he->next;
		nextv = he->sv;
		Loop* lp = loop;

		while (nextv != startv)
		{
			point[0] = nextv->coordinate[0] + d * dir[0];
			point[1] = nextv->coordinate[1] + d * dir[1];
			point[2] = nextv->coordinate[2] + d * dir[2];
			uphe = mev(nextv, point, lp);
			upv = uphe->ev;

			lp = mef(upprev, upv, loop, false);

			upprev = upv;
			he = he->next;
			nextv = he->sv;
		}
		mef(upprev, suphe->ev, loop, false);
		addLoopIntoFace(loop, loop->face->solid->faces);
	}
}//定义基本的欧拉操作
void addEdgeIntoSolid(Edge* edge, Solid*& solid)
{
	Edge* te = solid->edges;

	if (te == NULL)solid->edges = edge;
	else {
		while (te->next != NULL)te = te->next;
		te->next = edge;
		edge->pre = te;
	}
}
void addFaceIntoSolid(Face* face, Solid*& solid)
{
	Face* tface = solid->faces;
	if (tface == NULL)
	{
		solid->faces = face;
	}
	else
	{
		while (tface->next != NULL)tface = tface->next;
		tface->next = face;
		face->pre = tface;
	}
	face->solid = solid;

	face->id = solid->fnum;

	solid->fnum += 1;// increase the num of faces
}
void addLoopIntoFace(Loop* loop, Face* face)
{
	loop->face = face;

	
	if (face->out_lp == NULL)
	{
		face->out_lp = loop;
	}
	else if (face->out_lp == loop)
	{
		return;
	}
	else
	{
		Loop* tlp = face->inner_lp;
		if (tlp == NULL)
		{
			face->inner_lp = loop;
		}
		else
		{
			while (tlp->next != NULL)tlp = tlp->next;
			tlp->next = loop;
			loop->pre = tlp;
		}
		face->innum += 1;
	}
}
void construction()
{
	double Array[10][3] =
	{
		{0,0,0},
		{10,0,0},
		{10,20,0},
		{0,20,0},
		{1,1,0},
		{9,2,0},
		{2,9,0},
		{1,12,0},
		{9,12,0},
		{1,19,0},
	};
	Vertex* sv;
	Solid* solid = mvfs(Array[0], sv);
	HalfEdge* half1,*half2;
	half1 = mev(sv, Array[1], solid->faces->out_lp);
	half1 = mev(half1->ev, Array[2], solid->faces->out_lp);
	half1 = mev(half1->ev, Array[3], solid->faces->out_lp);
	Loop *loop1= mef(half1->ev, sv,solid->faces->out_lp, true);
	half1 = mev(sv, Array[4], solid->faces->out_lp);
	half2 = mev(half1->ev, Array[5], solid->faces->out_lp);
	half2 = mev(half2->ev, Array[6], solid->faces->out_lp);
	Loop *loop2 = mef(half2->ev,half1->ev,solid->faces->out_lp,true);
	Loop* loop3 = kemr(half1->sv, half2->next->ev, loop2);
	kfmrh(loop1->halfedges->brother->lp->face, loop2->halfedges->brother->lp->face);
	half1 = mev(sv, Array[7], solid->faces->out_lp);
	half2 = mev(half1->ev, Array[8], solid->faces->out_lp);
	half2 = mev(half2->ev, Array[9], solid->faces->out_lp);
	loop2 = mef(half2->ev, half1->ev, solid->faces->out_lp, true);
	loop2 = kemr(half1->sv, half2->next->ev, loop2);
	kfmrh(loop1->halfedges->brother->lp->face, loop2->halfedges->brother->lp->face);
	double direction[3] = { 0,0,1 };
	sweep(direction, 10);
	cout << solid->fnum << ' ' << solid->vnum << ' ' << solid->lnum << endl;
} //构建三维物体的半边数据结构
void construction1()
{
	double Array1[10][3];
	Vertex* sv;
	for (int k = 0; k < 10; k ++)
	{
		cout << "请输入底面第" << k+1  << "个点的坐标(注意按照逆时针方向依次输入外环（4个点）、内环（2*3个点）)" << endl;
		cin >> Array1[k][0]; 
		cin >> Array1[k][1];
		cin >> Array1[k][2];
	}
	Solid* solid = mvfs(Array1[0], sv);
	HalfEdge* half1, * half2;
	half1 = mev(sv, Array1[1], solid->faces->out_lp);
	half1 = mev(half1->ev, Array1[2], solid->faces->out_lp);
	half1 = mev(half1->ev, Array1[3], solid->faces->out_lp);
	Loop* loop1 = mef(half1->ev, sv, solid->faces->out_lp, true);
	half1 = mev(sv, Array1[4], solid->faces->out_lp);
	half2 = mev(half1->ev, Array1[5], solid->faces->out_lp);
	half2 = mev(half2->ev, Array1[6], solid->faces->out_lp);
	Loop* loop2 = mef(half2->ev, half1->ev, solid->faces->out_lp, true);
	Loop* loop3 = kemr(half1->sv, half2->next->ev, loop2);
	kfmrh(loop1->halfedges->brother->lp->face, loop2->halfedges->brother->lp->face);
	half1 = mev(sv, Array1[7], solid->faces->out_lp);
	half2 = mev(half1->ev, Array1[8], solid->faces->out_lp);
	half2 = mev(half2->ev, Array1[9], solid->faces->out_lp);
	loop2 = mef(half2->ev, half1->ev, solid->faces->out_lp, true);
	loop2 = kemr(half1->sv, half2->next->ev, loop2);
	kfmrh(loop1->halfedges->brother->lp->face, loop2->halfedges->brother->lp->face);
	double direction[3];
	int a;
	cout << "请输入扫掠的方向（单位向量）" << endl;
	cin >> direction[0] >> direction[1] >> direction[2];
	cout << "请输入在扫掠方向上的边长" <<endl;
	cin >> a;
	sweep(direction, a);
	cout << "构建的面的数量为：" << solid->fnum << endl;
	cout << "构建的点的数量为：" << solid->vnum << endl;
	cout<<"构建的环的数量为："<< solid->lnum << endl;
	cout << "三维模型生成中......" << endl;
} //构建三维物体的半边数据结构
void showloop()
{
	for (int i = 0; i < l_list.size(); i++)
	{
		cout << l_list[i]->id << endl;
		cout << l_list[i]->halfedges->sv->id << ' ';
		cout << l_list[i]->halfedges->next->sv->id << ' ';
		cout << l_list[i]->halfedges->next->next->sv->id << ' ';
		cout << l_list[i]->halfedges->next->next->next->sv->id << endl;
		cout << l_list[i]->face->innum;
		cout << endl << endl;
	}

}; 
void dispoint()
{
	for (int i = 0;i < v_list.size(); i++)
	{
		cout << v_list[i]->coordinate[0] << ' ' << v_list[i]->coordinate[1] << ' ' << v_list[i]->coordinate[2] << endl;
	}
}
void makeface(int a,int b,int c,int d)
{
	glVertex3f(v_list[a]->coordinate[0], v_list[a]->coordinate[1], v_list[a]->coordinate[2]);
	glVertex3f(v_list[b]->coordinate[0], v_list[b]->coordinate[1], v_list[b]->coordinate[2]);
	glVertex3f(v_list[c]->coordinate[0], v_list[c]->coordinate[1], v_list[c]->coordinate[2]);
	glVertex3f(v_list[d]->coordinate[0], v_list[d]->coordinate[1], v_list[d]->coordinate[2]);
}
void makeface1(int a, int b, int c)
{
	glVertex3f(v_list[a]->coordinate[0], v_list[a]->coordinate[1], v_list[a]->coordinate[2]);
	glVertex3f(v_list[b]->coordinate[0], v_list[b]->coordinate[1], v_list[b]->coordinate[2]);
	glVertex3f(v_list[c]->coordinate[0], v_list[c]->coordinate[1], v_list[c]->coordinate[2]);
}
void makeline(int a, int b)
{
	glVertex3f(v_list[a]->coordinate[0], v_list[a]->coordinate[1], v_list[a]->coordinate[2]);
	glVertex3f(v_list[b]->coordinate[0], v_list[b]->coordinate[1], v_list[b]->coordinate[2]);
}
void makeline3(Loop* loop)
{
	HalfEdge* Half = loop->halfedges;
	Vertex* sp = Half->sv;
	Vertex* ep = Half->ev;
	glBegin(GL_LINE_LOOP);
	glColor3f(1.0f, 1.0f, 1.0f);
	while (sp != ep)
	{
		glVertex3f(ep->coordinate[0], ep->coordinate[1], ep->coordinate[2]);
		Half = Half->next;
		ep = Half->ev;
	}
	glVertex3f(ep->coordinate[0], ep->coordinate[1], ep->coordinate[2]);
	glEnd();
}
void makeface3(Loop *loop)
{
	HalfEdge* Half = loop->halfedges;
	Vertex* sp = Half->sv;
	Vertex* ep = Half->ev;
	glBegin(GL_TRIANGLE_FAN);
	glColor3f(0.5f, 0.5f, 0.0f);
	while (sp != ep)
	{
		glVertex3f(ep->coordinate[0], ep->coordinate[1], ep->coordinate[2]);
		Half = Half->next;
		ep = Half->ev;
	}
	glVertex3f(ep->coordinate[0], ep->coordinate[1], ep->coordinate[2]);
	glEnd();
}
void makeSpecialFace1(Loop* loop)
{
	double midx1 = (v_list[10]->coordinate[0] + v_list[11]->coordinate[0]) / 2;
	double midy1 = (v_list[10]->coordinate[1] + v_list[11]->coordinate[1]) / 2;
	double midz1 = (v_list[10]->coordinate[2] + v_list[11]->coordinate[2]) / 2;
	double midx2 = (v_list[12]->coordinate[0] + v_list[13]->coordinate[0]) / 2;
	double midy2 = (v_list[12]->coordinate[1] + v_list[13]->coordinate[1]) / 2;
	double midz2 = (v_list[12]->coordinate[2] + v_list[13]->coordinate[2]) / 2;
	glBegin(GL_TRIANGLE_STRIP);
	glColor3f(0.4f, 0.4f, 0.0f);
	glVertex3f(midx1, midy1, midz1);
	glVertex3f(v_list[15]->coordinate[0], v_list[15]->coordinate[1], v_list[15]->coordinate[2]);
	glVertex3f(midx2, midy2, midz2);
	glVertex3f(v_list[16]->coordinate[0], v_list[16]->coordinate[1], v_list[16]->coordinate[2]);
	glVertex3f(v_list[12]->coordinate[0], v_list[12]->coordinate[1], v_list[12]->coordinate[2]);
	glVertex3f(v_list[14]->coordinate[0], v_list[14]->coordinate[1], v_list[14]->coordinate[2]);
	glVertex3f(v_list[11]->coordinate[0], v_list[11]->coordinate[1], v_list[11]->coordinate[2]);
	glVertex3f(v_list[15]->coordinate[0], v_list[15]->coordinate[1], v_list[15]->coordinate[2]);
	glVertex3f(midx1, midy1, midz1);
	glEnd();
	glBegin(GL_TRIANGLE_STRIP);
	glColor3f(0.4f, 0.4f, 0.0f);
	glVertex3f(midx1, midy1, midz1);
	glVertex3f(v_list[17]->coordinate[0], v_list[17]->coordinate[1], v_list[17]->coordinate[2]);
	glVertex3f(midx2, midy2, midz2);
	glVertex3f(v_list[19]->coordinate[0], v_list[19]->coordinate[1], v_list[19]->coordinate[2]);
	glVertex3f(v_list[13]->coordinate[0], v_list[13]->coordinate[1], v_list[13]->coordinate[2]);
	glVertex3f(v_list[18]->coordinate[0], v_list[18]->coordinate[1], v_list[18]->coordinate[2]);
	glVertex3f(v_list[10]->coordinate[0], v_list[10]->coordinate[1], v_list[10]->coordinate[2]);
	glVertex3f(v_list[17]->coordinate[0], v_list[17]->coordinate[1], v_list[17]->coordinate[2]);
	glVertex3f(midx1, midy1, midz1);
	glEnd();
}
void makeSpecialFace2(Loop* loop)
{
	double midx1 = (v_list[0]->coordinate[0] + v_list[3]->coordinate[0]) / 2;
	double midy1 = (v_list[0]->coordinate[1] + v_list[3]->coordinate[1]) / 2;
	double midz1 = (v_list[0]->coordinate[2] + v_list[3]->coordinate[2]) / 2;
	double midx2 = (v_list[2]->coordinate[0] + v_list[1]->coordinate[0]) / 2;
	double midy2 = (v_list[2]->coordinate[1] + v_list[1]->coordinate[1]) / 2;
	double midz2 = (v_list[2]->coordinate[2] + v_list[1]->coordinate[2]) / 2;
	glBegin(GL_TRIANGLE_STRIP);
	glColor3f(0.4f, 0.4f, 0.0f);
	glVertex3f(midx1, midy1, midz1);
	glVertex3f(v_list[6]->coordinate[0], v_list[6]->coordinate[1], v_list[6]->coordinate[2]);
	glVertex3f(v_list[0]->coordinate[0], v_list[0]->coordinate[1], v_list[0]->coordinate[2]);
	glVertex3f(v_list[4]->coordinate[0], v_list[4]->coordinate[1], v_list[4]->coordinate[2]);
	glVertex3f(v_list[1]->coordinate[0], v_list[1]->coordinate[1], v_list[1]->coordinate[2]);
	glVertex3f(v_list[5]->coordinate[0], v_list[5]->coordinate[1], v_list[5]->coordinate[2]);
	glVertex3f(midx2, midy2, midz2);
	glVertex3f(v_list[6]->coordinate[0], v_list[6]->coordinate[1], v_list[6]->coordinate[2]);
	glVertex3f(midx1, midy1, midz1);
	glEnd();
	glBegin(GL_TRIANGLE_STRIP);
	glColor3f(0.4f, 0.4f, 0.0f);
	glVertex3f(midx1, midy1, midz1);
	glVertex3f(v_list[7]->coordinate[0], v_list[7]->coordinate[1], v_list[7]->coordinate[2]);
	glVertex3f(midx2, midy2, midz2);
	glVertex3f(v_list[8]->coordinate[0], v_list[8]->coordinate[1], v_list[8]->coordinate[2]);
	glVertex3f(v_list[2]->coordinate[0], v_list[2]->coordinate[1], v_list[2]->coordinate[2]);
	glVertex3f(v_list[9]->coordinate[0], v_list[9]->coordinate[1], v_list[9]->coordinate[2]);
	glVertex3f(v_list[3]->coordinate[0], v_list[3]->coordinate[1], v_list[3]->coordinate[2]);
	glVertex3f(v_list[7]->coordinate[0], v_list[7]->coordinate[1], v_list[7]->coordinate[2]);
	glVertex3f(midx1, midy1, midz1);
	glEnd();

}
void display2()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);//清理屏幕颜色和深度
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1.0, 2.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	glLineWidth(2.0f);
	glColor3f(0.4f, 0.4f, 0.0f);
	makeSpecialFace1(l_list[0]);
	makeSpecialFace2(l_list[1]);
	for (int i = 0; i < l_list.size(); i++)
	{
	    if (l_list[i]->face->innum == 0)
		{
		    makeface3(l_list[i]);
		}
	    if (l_list[i]->face->innum == 0)
	    {
		    makeline3(l_list[i]);
	    }
	}
	glFlush();
}//显示
void display3()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);//清理屏幕颜色和深度
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	cout << "请输入你的观察视角的位置（例：5,4,3)" << endl;
	double a, b, c;
	cin >> a >> b >> c;
	gluLookAt(a, b, b, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	glLineWidth(2.0f);
	glColor3f(0.4f, 0.4f, 0.0f);
	makeSpecialFace1(l_list[0]);
	makeSpecialFace2(l_list[1]);
	for (int i = 0; i < l_list.size(); i++)
	{
		if (l_list[i]->face->innum == 0)
		{
			makeface3(l_list[i]);
		}
		if (l_list[i]->face->innum == 0)
		{
			makeline3(l_list[i]);
		}
	}
	glFlush();
}//显示
void shapetransform(int w, int h)
{
	glMatrixMode(GL_PROJECTION);       //设置为投影模式
	glLoadIdentity();
	glOrtho(-20.0, 20.0, -20.0, 20.0, -20.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}

int main(int argc, char** argv)
{
	cout << "默认or自己构造(1为默认，0自己构造)" << endl;
	bool a;
	cin >> a;
	if (a == 1)
	{
		construction();//默认设置构造
	}
	else
	{
		construction1();//自行构造
	}
	glutInit(&argc, argv); //来初始化GLUT库并同窗口系统对话协商
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(500, 500);
	glutInitDisplayMode(GLUT_RGBA);
	glutCreateWindow("3D widget");
	glutDisplayFunc(display2);//默认设置观看
	glutReshapeFunc(shapetransform);
	glutMainLoop();
	return 0;
}