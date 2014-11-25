/**********************************************************

                       GenFunctions

**********************************************************/


/*

    (c) Copyright 2002 by Michael J. Sanderson.

    Permission is granted to copy and use this program provided 
    no fee is charged for it and provided that this copyright 
    notice is not removed.
  

    This file is part of Pofad, (c) Copyright 2005 by Simon Joly.

    Pofad is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Pofad is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/


#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <ctype.h>

#define MAX_TOKEN_SIZE 10000		/* we've got room */

using namespace std;

void			doGenericAlert(char* errorMsg);
void			strtoupper(char *s);
void 			fatal(char *s);
