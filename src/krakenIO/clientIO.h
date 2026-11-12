/*---------------------------------*- C -*-----------------------------------*\
Copyright (C) 2026 VTT Technical Research Centre of Finland Ltd
-------------------------------------------------------------------------------
License
    This file is part of KrakenFOAM.

    KrakenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KrakenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with KrakenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Declarations of the C socket communication routines in clientIO.c.

\*---------------------------------------------------------------------------*/

#ifndef clientIO_H
#define clientIO_H

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

// Return values: 0 on success, -1 if sending or connecting failed, -2 if
// receiving failed and -3 if a received string does not fit in the buffer.
// The connection is closed if sending or receiving fails.

//- Connect to the server. A negative port number selects a unix socket
//  ./.socket_<-port>, otherwise a TCP connection to hostName:port is opened
int openConnection(const char* hostName, int64_t port);

//- Close the connection
int closeConnection(void);

//- Send length bytes from buffer
int sendBytes(const char* buffer, int64_t length);

//- Receive length bytes to buffer
int receiveBytes(char* buffer, int64_t length);

//- Send a 64-bit integer
int sendInt64(int64_t value);

//- Receive a 64-bit integer
int receiveInt64(int64_t* value);

//- Send a double
int sendDouble(double value);

//- Receive a double
int receiveDouble(double* value);

//- Send a string preceded by its length. The string is not null terminated.
int sendString(const char* str, int64_t length);

//- Receive a string preceded by its length to a buffer of size maxLength.
//  The string is not null terminated.
int receiveString(char* str, int64_t maxLength, int64_t* length);

//- Send an array of 64-bit integers
int sendInt64Array(const int64_t* array, int64_t length);

//- Receive an array of 64-bit integers
int receiveInt64Array(int64_t* array, int64_t length);

//- Send an array of doubles
int sendDoubleArray(const double* array, int64_t length);

//- Receive an array of doubles
int receiveDoubleArray(double* array, int64_t length);

#ifdef __cplusplus
}
#endif

#endif

/* ************************************************************************* */
