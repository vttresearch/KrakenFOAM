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
    C routines for socket communication with the Cerberus multiphysics
    driver. Used by Foam::krakenIO.

\*---------------------------------------------------------------------------*/

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32

#pragma comment(lib, "ws2_32.lib")
#define _WINSOCK_DEPRECATED_NO_WARNINGS
#include <winsock2.h>
#include <ws2tcpip.h>

#define SOCKET_LENGTH int
#define SEND_FLAGS 0
#define INTERRUPTED 0

#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/un.h>
#include <netdb.h>

#define SOCKET int
#define INVALID_SOCKET (-1)
#define SOCKET_ERROR (-1)
#define closesocket close
#define SOCKET_LENGTH size_t
#define INTERRUPTED (errno == EINTR)

// Prevent SIGPIPE from terminating the process if the server has closed the
// connection, the error is returned by send instead
#ifdef MSG_NOSIGNAL
    #define SEND_FLAGS MSG_NOSIGNAL
#else
    #define SEND_FLAGS 0
#endif

#endif

#include "clientIO.h"

#define PORT_STRING_SIZE 32

static SOCKET connection = INVALID_SOCKET;
static int isConnected = 0;


/***************************** PRIVATE FUNCTIONS *****************************/

static void printHostLookupError(const char* hostName)
{
    switch (h_errno)
    {
        case HOST_NOT_FOUND:
            fprintf(stderr, "Simulation host %s not found\n", hostName);
            break;
        case NO_ADDRESS:
            fprintf
            (
                stderr,
                "Invalid host address or no address found for simulation "
                "host %s\n",
                hostName
            );
            break;
        case NO_RECOVERY:
            fprintf(stderr, "A non-recoverable name server error occurred\n");
            break;
        case TRY_AGAIN:
            fprintf
            (
                stderr,
                "A temporary name server error occurred, host %s not found\n",
                hostName
            );
            break;
        default:
            fprintf
            (
                stderr,
                "An unknown error occurred while looking up the address of "
                "host %s\n",
                hostName
            );
    }
}


static void dropConnection(void)
{
    closesocket(connection);
    connection = INVALID_SOCKET;
    isConnected = 0;
}


/***************************** GLOBAL FUNCTIONS ******************************/

int openConnection(const char* hostName, int64_t port)
{
#ifdef _WIN32

    WSADATA wsaData;
    struct addrinfo* result = NULL;
    struct addrinfo* ptr = NULL;
    struct addrinfo hints;
    char portString[PORT_STRING_SIZE];
    int status;

    // TODO: Use hostName instead of the local host
    (void)hostName;

    if (port < 0)
    {
        fprintf(stderr, "Unix sockets on Windows not implemented yet.\n");
        return -1;
    }

    status = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (status != 0)
    {
        fprintf(stderr, "WSAStartup failed with error: %d\n", status);
        return -1;
    }

    ZeroMemory(&hints, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_protocol = IPPROTO_TCP;

    // Resolve the server address and port
    sprintf_s(portString, PORT_STRING_SIZE, "%" PRId64, port);
    status = getaddrinfo("127.0.0.1", portString, &hints, &result);
    if (status != 0)
    {
        fprintf(stderr, "getaddrinfo failed with error: %d\n", status);
        WSACleanup();
        return -1;
    }

    // Attempt to connect to an address until one succeeds
    connection = INVALID_SOCKET;
    for (ptr = result; ptr != NULL; ptr = ptr->ai_next)
    {
        connection =
            socket(ptr->ai_family, ptr->ai_socktype, ptr->ai_protocol);

        if (connection == INVALID_SOCKET)
        {
            fprintf
            (
                stderr,
                "socket failed with error: %d\n",
                WSAGetLastError()
            );
            freeaddrinfo(result);
            WSACleanup();
            return -1;
        }

        status = connect(connection, ptr->ai_addr, (int)ptr->ai_addrlen);
        if (status == SOCKET_ERROR)
        {
            fprintf
            (
                stderr,
                "connect failed with error: %d\n",
                WSAGetLastError()
            );
            closesocket(connection);
            connection = INVALID_SOCKET;
            continue;
        }
        break;
    }

    freeaddrinfo(result);

    if (connection == INVALID_SOCKET)
    {
        fprintf(stderr, "Unable to connect to server!\n");
        WSACleanup();
        return -1;
    }

#else

    struct sockaddr_storage address;
    socklen_t addressLength;
    int family;

    memset(&address, 0, sizeof(address));

    if (port < 0)
    {
        // Unix socket, no need to set up hostnames and IP addresses
        struct sockaddr_un* addressUnix = (struct sockaddr_un*)&address;

        family = AF_UNIX;
        addressUnix->sun_family = AF_UNIX;
        snprintf
        (
            addressUnix->sun_path,
            sizeof(addressUnix->sun_path),
            "./.socket_%" PRId64,
            -port
        );
        addressLength = sizeof(struct sockaddr_un);
    }
    else
    {
        struct hostent* host = gethostbyname(hostName);

        if (host == NULL)
        {
            printHostLookupError(hostName);
            return -1;
        }

        if (host->h_addrtype == AF_INET)
        {
            struct sockaddr_in* address4 = (struct sockaddr_in*)&address;

            family = AF_INET;
            address4->sin_family = AF_INET;
            address4->sin_port = htons((uint16_t)port);
            memcpy
            (
                &address4->sin_addr,
                host->h_addr_list[0],
                sizeof(struct in_addr)
            );
            addressLength = sizeof(struct sockaddr_in);
        }
        else if (host->h_addrtype == AF_INET6)
        {
            struct sockaddr_in6* address6 = (struct sockaddr_in6*)&address;

            family = AF_INET6;
            address6->sin6_family = AF_INET6;
            address6->sin6_port = htons((uint16_t)port);
            memcpy
            (
                &address6->sin6_addr,
                host->h_addr_list[0],
                sizeof(struct in6_addr)
            );
            addressLength = sizeof(struct sockaddr_in6);
        }
        else
        {
            fprintf(stderr, "Unknown address family for host %s\n", hostName);
            return -1;
        }
    }

    connection = socket(family, SOCK_STREAM, 0);
    if (connection == INVALID_SOCKET)
    {
        perror("Could not create socket");
        return -1;
    }

    if
    (
        connect(connection, (struct sockaddr*)&address, addressLength)
     == SOCKET_ERROR
    )
    {
        perror("Cannot connect to server");
        closesocket(connection);
        connection = INVALID_SOCKET;
        return -1;
    }

#endif

    isConnected = 1;
    return 0;
}


int closeConnection(void)
{
    if (connection != INVALID_SOCKET)
    {
        closesocket(connection);
        connection = INVALID_SOCKET;
    }

#ifdef _WIN32
    WSACleanup();
#endif

    isConnected = 0;
    return 0;
}


int sendBytes(const char* buffer, int64_t length)
{
    int64_t position = 0;

    if (!isConnected)
    {
        return -1;
    }

    while (position < length)
    {
        const int64_t count =
            send
            (
                connection,
                buffer + position,
                (SOCKET_LENGTH)(length - position),
                SEND_FLAGS
            );

        if (count < 0 && INTERRUPTED)
        {
            continue;
        }

        if (count <= 0)
        {
            dropConnection();
            return -1;
        }

        position += count;
    }

    return 0;
}


int receiveBytes(char* buffer, int64_t length)
{
    int64_t position = 0;

    if (!isConnected)
    {
        return -2;
    }

    while (position < length)
    {
        const int64_t count =
            recv
            (
                connection,
                buffer + position,
                (SOCKET_LENGTH)(length - position),
                0
            );

        if (count < 0 && INTERRUPTED)
        {
            continue;
        }

        if (count <= 0)
        {
            dropConnection();
            return -2;
        }

        position += count;
    }

    return 0;
}


int sendInt64(int64_t value)
{
    return sendBytes((const char*)&value, sizeof(int64_t));
}


int receiveInt64(int64_t* value)
{
    return receiveBytes((char*)value, sizeof(int64_t));
}


int sendDouble(double value)
{
    return sendBytes((const char*)&value, sizeof(double));
}


int receiveDouble(double* value)
{
    return receiveBytes((char*)value, sizeof(double));
}


int sendString(const char* str, int64_t length)
{
    const int status = sendInt64(length);

    if (status != 0)
    {
        return status;
    }

    return sendBytes(str, length*sizeof(char));
}


int receiveString(char* str, int64_t maxLength, int64_t* length)
{
    const int status = receiveInt64(length);

    if (status != 0)
    {
        return status;
    }

    if (*length < 0 || *length > maxLength)
    {
        fprintf
        (
            stderr,
            "Received string length %" PRId64 " exceeds the buffer size %"
            PRId64 "\n",
            *length,
            maxLength
        );
        return -3;
    }

    return receiveBytes(str, *length*sizeof(char));
}


int sendInt64Array(const int64_t* array, int64_t length)
{
    return sendBytes((const char*)array, length*sizeof(int64_t));
}


int receiveInt64Array(int64_t* array, int64_t length)
{
    return receiveBytes((char*)array, length*sizeof(int64_t));
}


int sendDoubleArray(const double* array, int64_t length)
{
    return sendBytes((const char*)array, length*sizeof(double));
}


int receiveDoubleArray(double* array, int64_t length)
{
    return receiveBytes((char*)array, length*sizeof(double));
}


/* ************************************************************************* */
