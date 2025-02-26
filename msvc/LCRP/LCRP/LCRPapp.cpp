// ParallaxBA.cpp : Defines the exported functions for the DLL application.

#include "stdafx.h"
#include "LCRP.h"
#include "../../../src/LCRP/LCRPImp.h"

LCRPapi ILCRP*	newLCRP()
{
	return (ILCRP*)(new CLCRP );
}

LCRPapi void    freeLCRP( ILCRP* ptr )
{
	free( ptr );
}