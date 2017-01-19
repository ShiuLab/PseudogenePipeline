
/* copyright (c) 1997 William R. Pearson */

/* used only in Mac versions for file selection */
/* should use navigation services if available */


#include <Dialogs.h>
#include <Fonts.h>
#include <Types.h>
#include <Gestalt.h>
#include <Resources.h>
#include <Controls.h>
#include <StandardFile.h>
#include <Files.h>
#include <Folders.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NIL nil
#define PauseID	301
#define ExitID	302
#define FileDID 204
#define SFileDID 205

void HiliteDlog(DialogPtr);

SFTypeList tlist={'TEXT',0L,0L,0L};

extern Point wpos;
	
FileDlog(prompt,freply)
	char *prompt;
	StandardFileReply *freply;
{
	Point dpos={-1,-1};
	if (GetResource('DLOG',SFileDID)==NIL) {
		fprintf(stderr," cannot load %d DLOG resource\n",SFileDID); exit(1);
		}
	CtoPstr(prompt);
	ParamText((StringPtr)prompt,"\p","\p","\p");
/*	SFPGetFile(wpos, (StringPtr)prompt, 0L,(short)1, tlist, 0L, freply, FileDID, NIL); */
	CustomGetFile(NIL,
				-1,
				nil,
				freply,
				SFileDID,
				dpos,
				nil,
				nil,nil,nil,nil);

	ParamText("\p","\p","\p","\p");
	PtoCstr((StringPtr)prompt);
	}
	
TFileDlog(prompt,freply,plist,nl)
	char *prompt;
	StandardFileReply *freply;
	SFTypeList plist;
	int nl;
{
	Point dpos={-1,-1};
	if (GetResource('DLOG',SFileDID)==NIL) {
		fprintf(stderr," cannot load %d TFile DLOG resource\n",SFileDID); exit(1);
		}
	CtoPstr(prompt);
	ParamText((StringPtr)prompt,"\p","\p","\p");
/*	SFPGetFile(wpos,(StringPtr)prompt,0L,(short)nl,plist,0L,freply,FileDID,NIL); */
	CustomGetFile(NIL,
				nl,
				plist,
				freply,
				SFileDID,
				dpos,
				nil,
				nil,nil,nil,nil);
	ParamText("\p","\p","\p","\p");
	PtoCstr((StringPtr)prompt);
	}

SFileDlog(prompt,freply)
	char *prompt;
	StandardFileReply *freply;
{
	Point dpos={-1,-1};

	if (GetResource('DLOG',SFileDID)==NIL) {
		fprintf(stderr," cannot load %d DLOG resource\n",SFileDID); exit(1);
		}

	CtoPstr(prompt);
	ParamText((StringPtr)prompt,"\p","\p","\p");

/* 	StandardGetFile(NIL,(short)1,tlist,freply); */
	CustomGetFile(NIL,
				-1,
				nil,
				freply,
				SFileDID,
				dpos,
				nil,
				nil,nil,nil,nil);
	ParamText("\p","\p","\p","\p");
	PtoCstr((StringPtr)prompt);
	}
	
STFileDlog(char *prompt, StandardFileReply *freply,
	   SFTypeList plist, int nl)
{
  Point dpos={-1,-1};

  if (GetResource('DLOG',SFileDID)==NIL) {
    fprintf(stderr," cannot load %d TFile DLOG resource\n",SFileDID); exit(1);
  }
  CtoPstr(prompt);
  ParamText((StringPtr)prompt,"\p","\p","\p");

  CustomGetFile(NIL,
		-1,
		nil,
		freply,
		SFileDID,
		dpos,
		nil,
		nil,nil,nil,nil);
  ParamText("\p","\p","\p","\p");
  PtoCstr((StringPtr)prompt);
}
	
PauseAlert(unsigned char *prompt)
{
  if (GetResource('DLOG',PauseID)==NIL) {
    fprintf(stderr," cannot load %d TFile DLOG resource\n",PauseID); exit(1);
  }
  CtoPstr((char *)prompt);
  ParamText(prompt,"\p","\p","\p");
  CautionAlert(PauseID,NULL);
  ParamText("\p","\p","\p","\p");
}

IntroDlog(int DlogID, unsigned char *prompt)
{
  short itemHit;
  DialogPtr DP;

  CtoPstr((char *)prompt);
  ParamText(prompt,"\p","\p","\p");

  if (GetResource('DLOG',DlogID)==NIL) {
    fprintf(stderr," cannot load %d Intro DLOG resource\n",DlogID); exit(1);
  }
  DP = GetNewDialog(DlogID,NULL,(WindowPtr)-1);
  ShowWindow(DP);
  SelectWindow(DP);
  HiliteDlog(DP);
	
  ModalDialog(0L,&itemHit);
  DisposeDialog(DP);
  ParamText("\p","\p","\p","\p");
  PtoCstr(prompt);
}

NIntroDlog(int DlogID,unsigned char *p0,unsigned char *p1,
	   unsigned char *p2,unsigned char *p3)
{
  short itemHit;
  DialogPtr DP;
  unsigned char *p;

  for (p=p0; *p; p++) if (*p=='\n') *p=' ';
  for (p=p1; *p; p++) if (*p=='\n') *p=' ';
  for (p=p2; *p; p++) if (*p=='\n') *p=' ';
  for (p=p2; *p; p++) if (*p=='\n') *p=' ';

  CtoPstr((char *)p0);
  CtoPstr((char *)p1);
  CtoPstr((char *)p2);
  CtoPstr((char *)p3);
  ParamText(p0,p1,p2,p3);

  if (GetResource('DLOG',DlogID)==NIL) {
    fprintf(stderr," cannot load %d Intro DLOG resource\n",DlogID); exit(1);
  }
  DP = GetNewDialog(DlogID,NULL,(WindowPtr)-1);
  ShowWindow(DP);
  SelectWindow(DP);
  HiliteDlog(DP);
	
  ModalDialog(0L,&itemHit);
  DisposeDialog(DP);
  ParamText("\p","\p","\p","\p");
  PtoCstr(p0);
  PtoCstr(p1);
  PtoCstr(p2);
  PtoCstr(p3);
}

void
HiliteDlog(DialogPtr DP)
{
  Rect tRect;
  short  tType;
  Handle tItem;

  SetPort(DP);
  GetDialogItem(DP,1,&tType,&tItem,&tRect);
  PenSize(3, 3);		/* Change pen to draw thick default outline */
  InsetRect(&tRect, -4, -4);   	/* Draw outside the button by 1 pixel */
  FrameRoundRect(&tRect, 16, 16); /* Draw the outline */
  PenSize(1, 1);		/* Restore the pen size to the default value */
}
