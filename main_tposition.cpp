// GVDB library
#include "gvdb.h"			
#include "gvdb_render.h"	// OpenGL rendering
using namespace nvdb;

// T Position utils
#include "main.h"			// window system 
#include "nv_gui.h"			// gui system
#include <GL/glew.h>

VolumeGVDB		gvdb;
CUmodule		cuCustom;
CUfunction		cuRaycastKernel;
int				gui_function;

enum SELECTION {
	NONE = 1,
	LOWERY = 2,
	UPPERY = 3,
	LOWERX = 4,
	UPPERX = 5,
	LOWERZ = 6,
	UPPERZ = 7,
	SHOW = 8
};

class Tposition : public NVPWindow {
public:
	virtual bool init();
	virtual void display();
	virtual void reshape(int w, int h);
	virtual void motion(int x, int y, int dx, int dy);
	virtual void keyboardchar(unsigned char key, int mods, int x, int y);
	virtual void mouse(NVPWindow::MouseButton button, NVPWindow::ButtonAction state, int mods, int x, int y);
	virtual void keyboard(KeyCode key, ButtonAction action, int mods, int x, int y);
	virtual void defineRays(Vector3DF origin, Vector3DF direction);
	virtual float3 getFloat3(Vector3DF origin);
	virtual void printVector3DF(Vector3DF vector);

	bool		LoadRAW(char* fname, Vector3DI res, int bpp);
	bool		ConvertToFloat(Vector3DI res, uchar* dat);
	void		Rebuild() { Rebuild(m_VolMax, m_dense); }
	void		Rebuild(Vector3DF vmax, bool dense);
	void		compute();

	void		start_guis(int w, int h);
	void		draw_topology();
	void		draw_elipsoid();

	// selection workflow
	void		draw_elipsoid_selection();
	void		draw_elipsoid_selectX(float xlower, float xUpper);
	void		draw_elipsoid_selectY(float ylower, float yUpper);
	void		draw_elipsoid_selectZ(float zlower, float zUpper);
	void		selectioWorkflow();
	void		selectionWorkflowTransition(int newState);

	Vector3DI	m_DataRes;
	int			m_DataBpp;			// 1=byte, 2=ushort, 4=float
	char*		m_DataBuf;
	Vector3DF	m_VolMax;

	// selection workflow
	Vector3DF	selectionLower;
	Vector3DF	selectionUpper;
	int			selection;
	float		volY;
	float		volX;
	float		volZ;

	int			m_gvdb_tex;
	int			mouse_down;
	int			screenWidth;
	int			screenHeight;
	bool		m_show_topo;
	bool		m_elipsoid;
	bool		m_dense;
	bool		m_function;
};
Tposition t_position_obj;

void handle_gui(int gui, float val)
{
	if (gui == gui_function) {
		t_position_obj.compute();
	}

	if (gui >= 1) {
		t_position_obj.Rebuild();
	}
	t_position_obj.postRedisplay();
}

bool Tposition::init()
{
	screenWidth = getWidth();
	screenHeight = getHeight();			// window width & height
	selectionLower = Vector3DF(0, 0, 0);
	selectionUpper = Vector3DF(0, 0, 0);
	selection = NONE;
	m_gvdb_tex = -1;
	mouse_down = -1;
	m_show_topo = false;
	m_elipsoid = false;
	m_dense = false;
	m_function = false;
	m_DataBuf = 0;

	init2D("arial");

	char* path = "../resource";

	// ------------- Initialize GVDB ---------------------------
	int devid = -1;
	// used to hold statics and other performance information
	gvdb.SetVerbose(true);
	// use the first valid devie found
	gvdb.SetCudaDevice(GVDB_DEV_FIRST);
	gvdb.Initialize();
	gvdb.AddPath(path);
	gvdb.AddPath(ASSET_PATH);
	gvdb.StartRasterGL();
	// ------------- Initialize GVDB ---------------------------

	// Load raw file into the CPU
	LoadRAW(path, Vector3DI(402, 402, 279), 1);		// This sets m_DataRes and m_DataBuf
	printf("Convert to float.\n");
	ConvertToFloat(m_DataRes, (uchar*)m_DataBuf); // Convert data to float

	// Transfer source data to GPU
	printf("Transfer data to GPU.\n");
	// load the data into the AUX
	gvdb.PrepareAux(AUX_DATA3D, m_DataRes.x*m_DataRes.y*m_DataRes.z, sizeof(float), false, false);
	DataPtr& aux3D = gvdb.getAux(AUX_DATA3D);
	// set the DATA from the CPU
	gvdb.SetDataCPU(aux3D, m_DataRes.x*m_DataRes.y*m_DataRes.z, m_DataBuf, 0, sizeof(float));
	// commit the data on the GPU
	gvdb.CommitData(aux3D);

	// Rebuild the data in GVDB	
	m_VolMax = Vector3DF(500, 500, 400);
	Rebuild(m_VolMax, m_dense);

	gvdb.getScene()->SetSteps(0.25f, 16.f, 0.25f);			// Set raycasting steps
	gvdb.getScene()->SetVolumeRange(0.25f, 0.0f, 1.0f);		// Set volume value range
	gvdb.getScene()->SetCutoff(0.005f, 0.005f, 0.f);
	gvdb.getScene()->SetShadowParams(0, 0, 0);
	gvdb.getScene()->LinearTransferFunc(0.0f, 0.8f, Vector4DF(0, 0, 0, 0), Vector4DF(1.0f, 0.627f, 0.478f, 0.5f));
	gvdb.getScene()->LinearTransferFunc(0.8f, 1.0f, Vector4DF(1.0f, 0.627f, 0.478f, 0.5f), Vector4DF(1.0f, 0.627f, 0.478f, 0.8f));
	gvdb.CommitTransferFunc();

	// Create Camera 
	Camera3D* cam = new Camera3D;
	cam->setFov(50.0);
	cam->setNearFar(.1, 5000);
	cam->setOrbit(Vector3DF(90, 20, 0), Vector3DF(128, 200, 128), 1800, 1.0);
	gvdb.getScene()->SetCamera(cam);

	// Create Light
	Light* lgt = new Light;
	lgt->setOrbit(Vector3DF(299, 57.3f, 0), Vector3DF(0, 0, 0), 1400, 1.0);
	gvdb.getScene()->SetLight(0, lgt);

	// Add render buffer
	nvprintf("Creating screen buffer. %d x %d\n", screenWidth, screenHeight);
	gvdb.AddRenderBuf(0, screenWidth, screenHeight, 4);

	// display the gui interaction options
	start_guis(screenWidth, screenHeight);
	return true;
}

void Tposition::compute() {
	// possible starting point for calculations triggered by user
}

void Tposition::start_guis(int w, int h)
{
	clearGuis();
	setview2D(w, h);
	guiSetCallback(handle_gui);
	// HINT:  Gui user interface options
	addGui(20, h - 30, 130, 20, "Topology", GUI_CHECK, GUI_BOOL, &m_show_topo, 0, 1);
	addGui(160, h - 30, 130, 20, "Elipsoid", GUI_CHECK, GUI_BOOL, &m_elipsoid, 0, 1);
	addGui(300, h - 30, 130, 20, "Dense", GUI_CHECK, GUI_BOOL, &m_dense, 0, 1);
	gui_function = addGui(440, h - 30, 130, 20, "Function", GUI_CHECK, GUI_BOOL, &m_function, 0, 100);
}

bool Tposition::ConvertToFloat(Vector3DI res, uchar* dat)
{
	float* vnew = (float*)malloc(res.x*res.y*res.z * sizeof(float));
	float* vdest = vnew;
	uchar* vsrc = dat;

	for (int n = 0; n < res.x*res.y*res.z; n++)
		*vdest++ = float(*vsrc++) / 256.0f;

	if (m_DataBuf != 0) free(m_DataBuf);
	m_DataBuf = (char*)vnew;
	m_DataRes = res;
	m_DataBpp = 4;

	return true;
}

bool Tposition::LoadRAW(char* fname, Vector3DI res, int bpp)
{
	// Load RAW in CPU memory
	char scnpath[1024];
	printf("Loading volume data.\n");

	if (!gvdb.getScene()->FindFile("out.raw", scnpath)) {
		printf("Cannot find pvm file %s.\n", fname);
		nverror();
	}
	printf("Loading RAW. %s\n", scnpath);

	char buf[1024];

	// Load RAW file into memory
	FILE* fp = fopen(scnpath, "rb");
	if (fp == 0) {
		nvprintf("Cannot read %s.\n", scnpath);
		return false;
	}
	m_DataRes = res;
	m_DataBpp = bpp;
	int sz = res.x*res.y*res.z*bpp;

	if (m_DataBuf != 0x0) free(m_DataBuf);
	m_DataBuf = (char*)malloc(sz);

	fread(m_DataBuf, 1, sz, fp);

	fclose(fp);
}

void Tposition::Rebuild(Vector3DF vmax, bool dense)
{
	gvdb.Clear();
	gvdb.DestroyChannels();

	// Configure VDB 
	printf("Configure GVDB.\n");
	gvdb.Configure(3, 3, 3, 2, 3);
	//gvdb.Configure(3, 3, 3, 2, 3);
	gvdb.SetVoxelSize(1, 1, 1);
	gvdb.SetChannelDefault(16, 16, 16);
	gvdb.AddChannel(0, T_FLOAT, 1);

	// Set volume transform
	Matrix4F xform;
	// NOTE: 
	// Transform represents the mapping from output space to input space (GVDB space to source data).
	// SRT = Scale, Rotate, Translate.  p' = S R T p
	xform.SRT(Vector3DF(1, 0, 0), Vector3DF(0, 1, 0), Vector3DF(0, 0, 1), Vector3DF(0, 0, 0), Vector3DF(1, 1, 1));

	// Activate volume
	printf("Activate GVDB volume.\n");
	Extents e;
	e = gvdb.ComputeExtents(1, Vector3DF(0, 0, 0), vmax);

	if (dense) {
		// Dense - Activate all bricks in the data volume
		// To load data densely, we simply active every brick in the volume extents.
		gvdb.ActivateRegion(0, e);
	}
	else {
		// Sparse - Only activate bricks above threshold value
		// To load data sparsely, we downsample the volume to get the average value at each brick.		
		// DownsampleCPU will downsample the input AUX volume and retrieve the values back to CPU.
		//gvdb.DownsampleCPU(xform, m_DataRes, AUX_DATA3D, e.ires, vmax, AUX_DOWNSAMPLED, Vector3DF(0, 1, 0), Vector3DF(0, 1, 0));
		gvdb.DownsampleCPU(xform, m_DataRes, AUX_DATA3D, e.ires, vmax, AUX_DOWNSAMPLED, Vector3DF(0, 1, 0), Vector3DF(0, 1, 0));

		// Then we only activate bricks whose average is above some threshold.
		gvdb.ActivateRegionFromAux(e, AUX_DOWNSAMPLED, T_FLOAT, 0.00001f);
	}

	gvdb.FinishTopology();
	gvdb.UpdateAtlas();
	gvdb.ClearAllChannels();

	// Resample data. 
	// The two vectors here represent the input and output value ranges.
	// During ConvertToFloat we already divided the source uchar by 256, so input data is already [0,1]. 
	// Output data, stored in GVDB, will also be [0,1]. Third value is reserved for future (ignored).
	printf("Resample.\n");
	gvdb.Resample(0, xform, m_DataRes, AUX_DATA3D, Vector3DF(0, 1, 0), Vector3DF(0, 1, 0));
}

void Tposition::reshape(int w, int h)
{
	// Resize the opengl screen texture
	glViewport(0, 0, w, h);
	createScreenQuadGL(&m_gvdb_tex, w, h);

	// Resize the GVDB render buffers
	gvdb.ResizeRenderBuf(0, w, h, 4);

	// Resize 2D UIs
	start_guis(w, h);

	postRedisplay();
}

void Tposition::draw_topology()
{
	Vector3DF clrs[10];
	clrs[0] = Vector3DF(0, 0, 1);			// blue
	clrs[1] = Vector3DF(0, 1, 0);			// green
	clrs[2] = Vector3DF(1, 0, 0);			// red
	clrs[3] = Vector3DF(1, 1, 0);			// yellow
	clrs[4] = Vector3DF(1, 0, 1);			// purple
	clrs[5] = Vector3DF(0, 1, 1);			// aqua
	clrs[6] = Vector3DF(1, 0.5, 0);		// orange
	clrs[7] = Vector3DF(0, 0.5, 1);		// green-blue
	clrs[8] = Vector3DF(0.7f, 0.7f, 0.7f);	// grey

	Camera3D* cam = gvdb.getScene()->getCamera();

	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	Vector3DF bmin, bmax;
	Node* node;
	for (int lev = 0; lev < 5; lev++) {				// draw all levels
		int node_cnt = gvdb.getNumNodes(lev);
		for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
			node = gvdb.getNodeAtLevel(n, lev);
			bmin = gvdb.getWorldMin(node);		// get node bounding box
			bmax = gvdb.getWorldMax(node);		// draw node as a box
			drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, clrs[lev].x, clrs[lev].y, clrs[lev].z, 1);
		}
	}

	end3D();										// end 3D drawing
}

void Tposition::draw_elipsoid() {
	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	int lev = 0;
	int node_cnt = gvdb.getNumNodes(lev);
	Vector3DF bmin, bmax;
	Node* node;
	//	for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
	for (int n = 0; n < 1; n++) {
		node = gvdb.getNodeAtLevel(n, lev);
		bmin = gvdb.getWorldMin(node);		// get node bounding box
		bmax = gvdb.getWorldMax(node);		// draw node as a box
		//drawLine3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 0, 1, 1);
		drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 0, 1, 1);
	}
	end3D();										// end 3D drawing
}

void Tposition::draw_elipsoid_selection() {
	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	int lev = 0;
	int node_cnt = gvdb.getNumNodes(lev);
	Vector3DF bmin, bmax;
	Node* node;
	//	for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
	for (int n = 0; n < node_cnt; n++) {
		node = gvdb.getNodeAtLevel(n, lev);
		bmin = gvdb.getWorldMin(node);		// get node bounding box
		bmax = gvdb.getWorldMax(node);		// draw node as a box
											//drawLine3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 0, 1, 1);
		if (bmin.x >= selectionLower.x && bmax.x <= selectionUpper.x &&
			bmin.y >= selectionLower.y && bmax.y <= selectionUpper.y &&
			bmin.z >= selectionLower.z && bmax.z <= selectionUpper.z) {
			drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 1, 0, 1);
		}
	}
	end3D();										// end 3D drawing
}

void Tposition::draw_elipsoid_selectX(float xLower, float xUpper) {
	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	int lev = 0;
	int node_cnt = gvdb.getNumNodes(lev);
	Vector3DF bmin, bmax;
	Node* node;
	//	for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
	for (int n = 0; n < node_cnt; n++) {
		node = gvdb.getNodeAtLevel(n, lev);
		bmin = gvdb.getWorldMin(node);		// get node bounding box
		bmax = gvdb.getWorldMax(node);		// draw node as a box
		boolean drawBoxStepOne = (xLower == xUpper && bmin.x <= xLower && bmax.x >= xUpper);
		boolean drawBoxStepTwo = (xLower != xUpper && bmin.x >= xLower && bmax.x <= xUpper);
		if (drawBoxStepOne || drawBoxStepTwo) {
			drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 1, 0, 1);
		}
	}
	end3D();										// end 3D drawing
}

void Tposition::draw_elipsoid_selectY(float yLower, float yUpper) {
	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	int lev = 0;
	int node_cnt = gvdb.getNumNodes(lev);
	Vector3DF bmin, bmax;
	Node* node;
	//	for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
	for (int n = 0; n < node_cnt; n++) {
		node = gvdb.getNodeAtLevel(n, lev);
		bmin = gvdb.getWorldMin(node);		// get node bounding box
		bmax = gvdb.getWorldMax(node);		// draw node as a box
		boolean drawBoxStepOne = (yLower == yUpper && bmin.y <=yLower && bmax.y >=yLower);
		boolean drawBoxStepTwo = (yLower != yUpper && bmin.y >= yLower && bmax.y <= yUpper);

		if (drawBoxStepOne || drawBoxStepTwo){
			drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 1, 0, 1);
		}
	}
	end3D();										// end 3D drawing
}

void Tposition::draw_elipsoid_selectZ(float zLower, float zUpper) {
	start3D(gvdb.getScene()->getCamera());		// start 3D drawing
	int lev = 0;
	int node_cnt = gvdb.getNumNodes(lev);
	Vector3DF bmin, bmax;
	Node* node;
	//	for (int n = 0; n < node_cnt; n++) {			// draw all nodes at this level
	for (int n = 0; n < node_cnt; n++) {
		node = gvdb.getNodeAtLevel(n, lev);
		bmin = gvdb.getWorldMin(node);		// get node bounding box
		bmax = gvdb.getWorldMax(node);		// draw node as a box
		boolean drawBoxStepOne = (zLower == zUpper && bmin.z <= zLower && bmax.z >= zUpper);
		boolean drawBoxStepTwo = (zLower != zUpper && bmin.z >= zLower && bmax.z <= zUpper);

		if (drawBoxStepOne || drawBoxStepTwo) {
			drawBox3D(bmin.x, bmin.y, bmin.z, bmax.x, bmax.y, bmax.z, 0, 1, 0, 1);
		}
	}
	end3D();										// end 3D drawing
}

void Tposition::display()
{
	int w = getWidth(), h = getHeight();			// window width & height

	clearScreenGL();

	// Render volume
	gvdb.TimerStart();
	gvdb.Render(SHADE_TRILINEAR, 0, 0);    // last value indicates render buffer for depth input
	//gvdb.Render( SHADE_VOLUME, 0, 0 );    // last value indicates render buffer for depth input
	//float rtime = gvdb.TimerStop();
	//nvprintf("Render volume. %6.3f ms\n", rtime);

	// Copy GVDB output buffer into OpenGL texture	
	gvdb.ReadRenderTexGL(0, m_gvdb_tex);

	// Render to screen
	renderScreenQuadGL(m_gvdb_tex);

	if (m_show_topo) draw_topology();			// Draw GVDB topology 

	if (m_elipsoid) {
		volY = (1.0f - ((float)getCurY() / screenHeight)) * m_VolMax.y;
		volX = ((float)getCurX() / screenWidth) * m_VolMax.x;
		volZ = (1.0f - ((float)getCurX() / screenWidth)) * m_VolMax.x;
		this->selectioWorkflow();
	}

	draw3D();									// Render the 3D drawing groups
	drawGui(0);
	draw2D();
}

void Tposition::selectioWorkflow() {
	switch (selection) {
	case SELECTION::LOWERY: {
		draw_elipsoid_selectY(volY, volY);
	}break;
	case SELECTION::UPPERY: {
		draw_elipsoid_selectY(selectionLower.y, volY);
	}break;
	case SELECTION::LOWERX: {
		draw_elipsoid_selectX(volX, volX);
	}break;
	case SELECTION::UPPERX: {
		draw_elipsoid_selectX(selectionLower.x, volX);
	}break;
	case SELECTION::LOWERZ: {
		draw_elipsoid_selectZ(volZ, volZ);
	}break;
	case SELECTION::UPPERZ: {
		draw_elipsoid_selectZ(selectionLower.z, volZ);
	}break;
	case SELECTION::SHOW: {
		draw_elipsoid_selection();
	}break;
	}
}

void Tposition::selectionWorkflowTransition(int newState) {
	switch (newState) {
	case SELECTION::UPPERY: {
		selectionLower.y = volY-5;
	}break;
	case SELECTION::LOWERX: {
		selectionUpper.y = volY;
	}break;
	case SELECTION::UPPERX: {
		selectionLower.x = volX-5;
	}break;
	case SELECTION::LOWERZ: {
		selectionUpper.x = volX;
	}break;
	case SELECTION::UPPERZ: {
		selectionLower.z = volZ - 5;
	}break;
	case SELECTION::SHOW: {
		selectionUpper.z = volZ;
	}break;
	}
}

void Tposition::keyboardchar(unsigned char key, int mods, int x, int y)
{
	switch (key) {
	case '1':	m_show_topo = !m_show_topo;	break;
	};
}

void Tposition::motion(int x, int y, int dx, int dy)
{
	// Get camera for GVDB Scene
	Camera3D* cam = gvdb.getScene()->getCamera();
	Light* lgt = gvdb.getScene()->getLight();
	bool shift = (getMods() & NVPWindow::KMOD_SHIFT);		// Shift-key to modify light

	postRedisplay();

	switch (mouse_down) {
	case NVPWindow::MOUSE_BUTTON_LEFT: {
		// Adjust orbit angles
		Vector3DF angs = (shift ? lgt->getAng() : cam->getAng());
		angs.x += dx * 0.2f;
		angs.y -= dy * 0.2f;
		if (shift)	lgt->setOrbit(angs, lgt->getToPos(), lgt->getOrbitDist(), lgt->getDolly());
		else		cam->setOrbit(angs, cam->getToPos(), cam->getOrbitDist(), cam->getDolly());
		postRedisplay();
	} break;

	case NVPWindow::MOUSE_BUTTON_MIDDLE: {
		// Adjust target pos		
		cam->moveRelative(float(dx) * cam->getOrbitDist() / 1000, float(-dy) * cam->getOrbitDist() / 1000, 0);
		postRedisplay();
	} break;

	case NVPWindow::MOUSE_BUTTON_RIGHT: {
		// Adjust dist
		float dist = (shift ? lgt->getOrbitDist() : cam->getOrbitDist());
		dist -= dy;
		if (shift)	lgt->setOrbit(lgt->getAng(), lgt->getToPos(), dist, cam->getDolly());
		else		cam->setOrbit(cam->getAng(), cam->getToPos(), dist, cam->getDolly());
		postRedisplay();
	} break;
	}
}

void Tposition::mouse(NVPWindow::MouseButton button, NVPWindow::ButtonAction state, int mods, int x, int y)
{
	if (guiHandler(button, state, x, y)) return;

	// Track when we are in a mouse drag
	mouse_down = (state == NVPWindow::BUTTON_PRESS) ? button : -1;
}

void Tposition::keyboard(KeyCode key, ButtonAction action, int mods, int x, int y) {
	if ((action == BUTTON_PRESS && (key == KEY_LEFT_CONTROL || key == KEY_RIGHT_CONTROL)) && m_elipsoid) {
		selection++;
		if (selection == 9) {
			selection = 0;
		}
		nvprintf("Control pressed: selection %d\n", selection);
		selectionWorkflowTransition(selection);
	}
}

void Tposition::defineRays(Vector3DF origin, Vector3DF direction) {
	struct ALIGN(16) ScnRay {
		float3 hit; // hit point
		float3 normal; // hit normal
		float3 orig; // ray origin
		float3 dir; // ray direction
		uint clr; // ray color
		uint pnode;// internal
		uint pndx; // internal
	};
	float3 screenpositon = getFloat3(origin);
	float3 rayDirection = getFloat3(direction);
	DataPtr m_rays;
	int numRays = 1;
	gvdb.AllocData(m_rays, numRays, sizeof(ScnRay));
	// get CPU pointer to first ray
	ScnRay* ray = (ScnRay*)gvdb.getDataPtr(0, m_rays);
	for (int n = 0; n < numRays; n++) {
		ray->orig = screenpositon;
		ray->dir = rayDirection;
		ray++; // next ray
	}
	gvdb.CommitData(m_rays);

	// Raytrace
	gvdb.Raytrace(m_rays, 0, SHADE_TRILINEAR, 0, 0);

	gvdb.RetrieveData(m_rays);

	ScnRay* rayRetrieve = (ScnRay*)gvdb.getDataPtr(0, m_rays);
	for (int n = 0; n < numRays; n++) {
		// read data of the ray
		nvprintf("Hit at: x: %f, y: %f, z: %f\n", rayRetrieve->hit.x, rayRetrieve->hit.y, rayRetrieve->hit.z);
	}
}

float3 Tposition::getFloat3(Vector3DF origin) {
	float3 result;
	result.x = origin.x;
	result.y = origin.y;
	result.z = origin.z;
	return result;
}

void Tposition::printVector3DF(Vector3DF vector) {
	nvprintf("X: %f, Y: %f, Z: %f\n", vector.x, vector.y, vector.z);
}

int sample_main(int argc, const char** argv)
{
	return t_position_obj.run("Transform to T-position", "t-position", argc, argv, 1280, 720, 4, 5);
}

void sample_print(int argc, char const *argv)
{
}
