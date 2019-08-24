/*
** (c) 1996-2000 The Regents of the University of California (through
** E.O. Lawrence Berkeley National Laboratory), subject to approval by
** the U.S. Department of Energy.  Your use of this software is under
** license -- the license agreement is attached and included in the
** directory as license.txt or you may contact Berkeley Lab's Technology
** Transfer Department at TTD@lbl.gov.  NOTICE OF U.S. GOVERNMENT RIGHTS.
** The Software was developed under funding from the U.S. Government
** which consequently retains certain rights as follows: the
** U.S. Government has been granted for itself and others acting on its
** behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, and perform publicly
** and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of
** Energy, and subject to any subsequent five (5) year renewals, the
** U.S. Government is granted for itself and others acting on its behalf
** a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, distribute copies to
** the public, perform publicly and display publicly, and to permit
** others to do so.
*/

//
// $Id: BLThread.cpp,v 1.32 2002/03/26 22:14:33 car Exp $
//

#include <winstd.H>
#include <BoxLib.H>
#include <Thread.H>

#ifdef WIN32
#define _WIN32_WINNT 0x0400
#include <windows.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

#include <iostream>
#include <limits>

#if defined(BL_OLD_STL)
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <stdlib.h>
#else
#include <cstdio>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#endif

//#if defined(BL_OSF1)
//extern "C" int usleep (useconds_t);
//#endif

namespace
{
    const char*
    the_message_string(const char* file, int line, const char* call, int status = 0)
    {
	//
	// Should be large enough.
	//
	const int DIM = 1024;
	static char buf[DIM];
	if ( status )
	{
	    std::sprintf(buf, "BoxLib Thread Error: File %s, line %d, %s: %d",
			 file, line, call, status);
	}
	else
	{
	    std::sprintf(buf, "BoxLib Thread Error: File %s, line %d, %s",
			 file, line, call);
	}
	buf[DIM-1] = '\0';		// Just to be safe.
	return buf;
    }

}

namespace BoxLib
{
    void
    Thread_Error(const char* file, int line, const char* call, int status = 0)
    {
	Error(the_message_string(file, line, call,status));
    }

}

#define THREAD_REQUIRE(x)						\
do									\
{									\
  if ( int status = (x) )						\
    {									\
      BoxLib::Thread_Error(__FILE__, __LINE__, #x, status); 		\
    }									\
}									\
while ( false )

#define THREAD_ASSERT(x)						\
do									\
{									\
  if ( !(x) )								\
    {									\
      BoxLib::Thread_Error(__FILE__,__LINE__,#x );			\
    }									\
}									\
while ( false )

void
Thread::sleep(const BoxLib::Time& spec_)
{
#ifdef WIN32
#else
   ::sleep(spec_.as_long());
#endif
}


//
// Barrier
//

Barrier::Barrier(int i)
    : count(0), n_sleepers(0), releasing(false)
{
    init(i);
}

void
Barrier::init(int i)
{
    THREAD_ASSERT( !releasing );
    THREAD_ASSERT( n_sleepers == 0 );
    count = i;
}

void
Barrier::wait()
{
    bool release = false;
    lock();
    // If previous cycle still releasing, wait
    // THREAD_ASSERT ( !releasing );
    while ( releasing )
    {
	ConditionVariable::wait();
    }
    if ( ++n_sleepers == count )
    {
	release = releasing = true;
    }
    else
    {
	// A poor thread cancelation Site
	Thread::CancelState tmp = Thread::setCancelState(Thread::Disable);
	while ( !releasing )
	{
	    ConditionVariable::wait();
	}
	Thread::setCancelState(tmp);
    }
    if ( --n_sleepers == 0 )
    {
	releasing = false;
	release = true;             // Wake up waiters (if any) for next cycle
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}


//
// Semaphore
//

Semaphore::Semaphore(int val_)
    : value(val_)
{
}

void
Semaphore::wait()
{
    lock();
    while ( value == 0 )
    {
	ConditionVariable::wait();
    }
    value--;
    unlock();
}

bool
Semaphore::trywait()
{
    lock();
    if ( value == 0 )
    {
	unlock();
	return false;
    }
    value--;
    unlock();
    return true;
}

void
Semaphore::post()
{
    lock();
    value++;
    unlock();
    signal();
}


//
//
//

SemaphoreB::SemaphoreB(int val_)
    : val(val_)
{}

int
SemaphoreB::down()
{
    lock();
    while (val <= 0)
    {
	wait();
    }
    int t = --val;
    unlock();

    return t;
}

int
SemaphoreB::up()
{
    lock();
    int t = ++val;
    unlock();
    signal();
    return t;
}

int
SemaphoreB::decrement()
{
    lock();
    int t = --val;
    unlock();
    return t;
}

int
SemaphoreB::value()
{
    lock();
    int t = val;
    unlock();
    return t;
}

//
// SingleBarrier
//

SingleBarrier::SingleBarrier(int i)
    : count(i), n_posters(0), n_waiters(0), releasing(false)
{
}

void
SingleBarrier::wait()
{
    bool release = false;
    lock();
    n_waiters++;
    while ( !releasing )
    {
	ConditionVariable::wait();
    }
    if ( --n_waiters == 0 )
    {
	releasing = false;
	release = true;             // Wake up waiters (if any) for next cycle
	n_posters=0;
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}

void
SingleBarrier::post()
{
    bool release = false;
    lock();
    // If previous cycle still releasing, wait
    while ( releasing )
    {
	ConditionVariable::wait();
    }
    if ( ++n_posters == count )
    {
	releasing = true;
	release = true;             // Wake up waiters (if any) for next cycle
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}


//
//Gate
//

Gate::Gate()
    : closed(true)
{
}

void
Gate::open()
{
    lock();
    closed = false;
    broadcast();
    unlock();
}

void
Gate::close()
{
    lock();
    closed = true;
    unlock();
}

void
Gate::release()
{
    broadcast();
}

void
Gate::wait()
{
    lock();
    while ( closed )
    {
	ConditionVariable::wait();
    }
    unlock();
}


//
// Lock<Semaphore> specialization
//

Lock<Semaphore>::Lock(Semaphore& sem_)
    : sem(sem_)
{
    sem.wait();
}

Lock<Semaphore>::~Lock()
{
    sem.post();
}

void
Thread::exit(void*)
{
    std::exit(0);
}

Thread::CancelState
Thread::setCancelState(CancelState)
{
    return Enable;
}


//
//
//

class ThreadSpecificData<void>::Implementation
{
public:
    Implementation(void (*tsd_)(void*));
    ~Implementation();
    void* set(const void* v_);
    void* get() const;
    void* v;
    void (*tsd)(void*);
};

ThreadSpecificData<void>::Implementation::Implementation(void (*tsd_)(void*))
    : v(0), tsd(tsd_)
{
}

ThreadSpecificData<void>::Implementation::~Implementation()
{
}

void*
ThreadSpecificData<void>::Implementation::set(const void* v_)
{
    return v = const_cast<void*>(v_);
}

void*
ThreadSpecificData<void>::Implementation::get() const
{
    return v;
}


ThreadSpecificData<void>::ThreadSpecificData(void (*tsd_)(void*))
{
  m_impl = new Implementation(tsd_);
}

ThreadSpecificData<void>::~ThreadSpecificData()
{
  delete m_impl;
}

void*
ThreadSpecificData<void>::set(const void* v_)
{
    return m_impl->set(v_);
}

void*
ThreadSpecificData<void>::get() const
{
    return m_impl->get();
}

// FuctioinThread
FunctionThread::FunctionThread(Thread_Function func, void* arg_, DetachState st, int stacksize)
{
    func(arg_);
}

FunctionThread::~FunctionThread()
{
    detach();
}

void*
FunctionThread::join() const
{
    return 0;
}

void
FunctionThread::detach() const
{
}

