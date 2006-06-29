// C++ informative line for the emacs editor: -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
 * access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _H5Exception_H
#define _H5Exception_H

#include <string>

#ifndef H5_NO_NAMESPACE
namespace H5 {
#ifndef H5_NO_STD
    using namespace std;
#endif  // H5_NO_STD
#endif

class H5_DLLCPP Exception {
   public:
	// Creates an exception with a function name where the failure occurs
	// and an optional detailed message
	Exception(const string func_name, const string message = DEFAULT_MSG);

	// Returns a character string that describes the error specified by
	// a major error number.
	string getMajorString( H5E_major_t major_num ) const;

	// Returns a character string that describes the error specified by
	// a minor error number.
	string getMinorString( H5E_minor_t minor_num ) const;

	// Returns the detailed message set at the time the exception is thrown
	string getDetailMsg() const;
	const char* getCDetailMsg() const;	// C string of detailed message
	string getFuncName() const;	// function name as a string object
	const char* getCFuncName() const;	// function name as a char string

	// Turns on the automatic error printing.
	static void setAutoPrint( H5E_auto_t& func, void* client_data);

	// Turns off the automatic error printing.
	static void dontPrint();

	// Retrieves the current settings for the automatic error stack
	// traversal function and its data.
	static void getAutoPrint( H5E_auto_t& func, void** client_data);

	// Clears the error stack for the current thread.
	static void clearErrorStack();

	// Walks the error stack for the current thread, calling the
	// specified function.
	static void walkErrorStack( H5E_direction_t direction,
				H5E_walk_t func, void* client_data);

	// Prints the error stack in a default manner.
	virtual void printError( FILE* stream = NULL ) const;

	// Default constructor
	Exception();

	// copy constructor
	Exception( const Exception& orig);

	// virtual Destructor
	virtual ~Exception();

   private:
// Because 'string' is not instantiated at compilation time, this
// warning is displayed when building DLL; but the class is exported
// so the warning is harmless
#if defined(WIN32)
#pragma warning(disable: 4251)
#endif
	string detail_message;
	string func_name;

   protected:
        // Default value for detail_message
        static const string DEFAULT_MSG;
};

class H5_DLLCPP FileIException : public Exception {
   public:
	FileIException( const string func_name, const string message = DEFAULT_MSG);
	FileIException();
	virtual ~FileIException();
};

class H5_DLLCPP GroupIException : public Exception {
   public:
	GroupIException( const string func_name, const string message = DEFAULT_MSG);
	GroupIException();
	virtual ~GroupIException();
};

class H5_DLLCPP DataSpaceIException : public Exception {
   public:
	DataSpaceIException(const string func_name, const string message = DEFAULT_MSG);
	DataSpaceIException();
	virtual ~DataSpaceIException();
};

class H5_DLLCPP DataTypeIException : public Exception {
   public:
	DataTypeIException(const string func_name, const string message = DEFAULT_MSG);
	DataTypeIException();
	virtual ~DataTypeIException();
};

class H5_DLLCPP PropListIException : public Exception {
   public:
	PropListIException(const string func_name, const string message = DEFAULT_MSG);
	PropListIException();
	virtual ~PropListIException();
};

class H5_DLLCPP DataSetIException : public Exception {
   public:
	DataSetIException(const string func_name, const string message = DEFAULT_MSG);
	DataSetIException();
	virtual ~DataSetIException();
};

class H5_DLLCPP AttributeIException : public Exception {
   public:
	AttributeIException(const string func_name, const string message = DEFAULT_MSG);
	AttributeIException();
	virtual ~AttributeIException();
};

class H5_DLLCPP ReferenceException : public Exception {
   public:
	ReferenceException(const string func_name, const string message = DEFAULT_MSG);
	ReferenceException();
	virtual ~ReferenceException();
};

class H5_DLLCPP LibraryIException : public Exception {
   public:
	LibraryIException(const string func_name, const string message = DEFAULT_MSG);
	LibraryIException();
	virtual ~LibraryIException();
};

class H5_DLLCPP IdComponentException : public Exception {
   public:
	IdComponentException(const string func_name, const string message = DEFAULT_MSG);
	IdComponentException();
	virtual ~IdComponentException();
};

#ifndef H5_NO_NAMESPACE
}
#endif

#endif // _H5Exception_H
