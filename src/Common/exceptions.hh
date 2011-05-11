//  A. L. Delcher
//
//  File:  exceptions.hh
//
//  Last Modified:  13 June 2005
//
//  Include file for exception types



#ifndef  __EXCEPTIONS_HH_INCLUDED
#define  __EXCEPTIONS_HH_INCLUDED


// Stolen from AMOS exceptions code



const std :: string  NULL_STRING = "";  //!< null string



//================================================ Exception_t =================
//! \brief The base exception class
//!
//! All other exceptions will be derived from this class, so catching for
//! this class should effectively catch all exceptions.
//!
//==============================================================================

class Exception_t  :  public std :: exception
{

private:

  std :: string  what_m;    //!< description of exception
  int  line_m;              //!< line number of exception
  std ::  string file_m;    //!< file name of exception


public:

  //---------------------------------------------- Exception_t -----------------
  //! \brief Informative constructor
  //!
  //! \param what Brief description of the exception
  //! \param line Line number of the exception
  //! \param file File name of the exception
  //!
  Exception_t (const std :: string & what,
	       int line = 0,
	       const std :: string & file = NULL_STRING)
    : what_m (what), line_m (line), file_m (file)
  { }


  //---------------------------------------------- ~Exception_t ----------------
  //! \brief Default destructor
  //!
  ~Exception_t ( ) throw()
  { }


  //---------------------------------------------- file ------------------------
  //! \brief Returns the file (if available) of the exception
  //!
  virtual const char * file ( ) const
  {
    return file_m . c_str( );
  }


  //---------------------------------------------- line ------------------------
  //! \brief Returns the line number (if available) of the exception
  //!
  virtual int line ( ) const
  {
    return line_m;
  }


  //---------------------------------------------- what ------------------------
  //! \brief Returns a short description (if available) of the exception
  //!
  virtual const char * what ( ) const throw( )
  {
    return what_m . c_str( );
  }

};



//----------------------------------------------------- operator<< -------------
//! \brief Dump Exception_t info to an ostream
//!
inline std :: ostream & operator<< (std :: ostream & out, const Exception_t & e)
{
  out << "WHAT: " << e . what( ) << std::endl;
  out << "LINE: " << e . line( ) << std::endl;
  out << "FILE: " << e . file( ) << std::endl;
  return out;
}


//----------------------------------------------------- operator<< -------------
//! \brief Dump exception info to an ostream
//!
inline std::ostream & operator<< (std::ostream & out, const std::exception & e)
{
  out << "WHAT: " << e . what( ) << std::endl;
  return out;
}



//-- Helpful exception throw macros
#define SIMPLE_THROW(A) throw Exception_t (A, __LINE__, __FILE__)


#endif // #ifndef __EXCEPTIONS_HH_INCLUDED
