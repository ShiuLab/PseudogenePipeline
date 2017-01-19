
/* Copyright 1995 William R. Pearson */

/* used only in Mac versions to provide mac multitasking */

#include <stdlib.h>

#ifdef __MWERKS__
#include <sioux.h>
#endif

#define SLEEP			2L
#define NIL_MOUSE_REGION	0L

#define WNE_TRAP_NUM	0x60
#define UNIMPL_TRAP_NUM	0x9F
#define SUSPEND_RESUME_BIT	0x0001
#define ACTIVATING		1
#define RESUMING		1

Boolean		gDone, gWNEImplemented=0;
EventRecord	gTheEvent;
Rect		gDragRect, gSizeRect;

void
InitEvent()
{
	gWNEImplemented=(NGetTrapAddress(WNE_TRAP_NUM,ToolTrap)!=
		NGetTrapAddress(UNIMPL_TRAP_NUM,ToolTrap));
	}


#define hiword(x)		(((short *) &(x))[0])
#define loword(x)		(((short *) &(x))[1])
static MenuHandle aMenu;

/*
ChkEvent()
{}
*/

#ifdef TPLOT
extern WindowPtr gDrawWindow;
extern PicHandle aPic;
#endif

static long checkTime=0;

void
ChkEvent()
{
	EventRecord event;
	WindowPeek wp;
	Boolean gotEvent, SIOUXDidEvent;
	long choice;
	Str255 buf;
	
	if (TickCount() < checkTime) return;
	checkTime = TickCount()+60L;

	if (gWNEImplemented)
		gotEvent=WaitNextEvent(everyEvent-diskMask,&event,SLEEP,NIL_MOUSE_REGION);
	else {
		SystemTask();
		gotEvent=GetNextEvent(everyEvent-diskMask,&event);
		}

	if (gotEvent) SIOUXDidEvent=SIOUXHandleOneEvent(&event);
	if (SIOUXDidEvent) return;

	if (event.what == nullEvent) {
		if (FrontWindow() == 0) InitCursor();
		return;
		}
	
	if (SystemEvent(&event)) return;

	if (event.what == mouseDown) {
		switch (FindWindow(event.where, (WindowPtr *)&wp)) {
			case inMenuBar:
				InitCursor();
				choice = MenuSelect(event.where);
				goto doMenu;
			case inDrag :
				DragWindow((WindowPtr)wp, event.where, &gDragRect);
				break;
			case inSysWindow:
				SystemClick(&event, (WindowPtr)wp);
				break;
			}
		}
	
	return;

doMenu:	
	switch (hiword(choice)) {
		case 1:
			GetMenuItemText(aMenu, loword(choice), buf);
			OpenDeskAcc(buf);
			break;
		case 2:
			exit(0);
		
		case 3:	
			SystemEdit(loword(choice) - 1);
			break;
	}
	HiliteMenu(0);
}

#ifdef TPLOT

Waitkey(keyval)
	int keyval;
{
	int key;
	EventRecord event;
	WindowPeek wp;
	long choice;
	Str255 buf;
	
	SystemTask();
	if (gWNEImplemented)
		WaitNextEvent(everyEvent-diskMask,&event,SLEEP,NIL_MOUSE_REGION);
	else {
		SystemTask();
		GetNextEvent(everyEvent-diskMask,&event);
		}

			
	InitCursor();
	if (event.what == nullEvent) {
		return 0;
		}
	
	if (SystemEvent(&event)) return 0;

	if (event.what == updateEvt) {
		if ((WindowPtr)event.message == gDrawWindow) {
			BeginUpdate((WindowPtr)event.message);
			DrawPicture(aPic,&gDrawWindow->portRect);
			EndUpdate((WindowPtr)event.message);
			}
		else {
			BeginUpdate((WindowPtr)event.message);
			EndUpdate((WindowPtr)event.message);
			}
		return 0;
		}

	if (event.what == keyDown) return 1;
	if (event.what == mouseDown) {
		switch (FindWindow(event.where, (WindowPtr *)&wp)) {
			case inMenuBar:
				InitCursor();
				choice = MenuSelect(event.where);
				goto doMenu;
			case inDrag :
				DragWindow((WindowPtr)wp, event.where, &gDragRect);
				break;
			case inSysWindow:
				SystemClick(&event, (WindowPtr)wp);
				break;
			case inGoAway :
				return 1;
			case inContent:
				SelectWindow((WindowPtr)wp);
				SetPort(gDrawWindow);
				DrawPicture(aPic,&gDrawWindow->portRect);
				break;
			}
		}
	
	return 0;

doMenu:	
	switch (hiword(choice)) {
		case 1:
			GetItem(aMenu, loword(choice), buf);
			OpenDeskAcc(buf);
			break;
		case 2:
			return 1;
		
		case 3:	
			SystemEdit(loword(choice) - 1);
			break;
	}
	HiliteMenu(0);
	return 0;
}
#endif

			
