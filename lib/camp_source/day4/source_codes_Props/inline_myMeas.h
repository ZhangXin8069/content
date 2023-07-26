// -*- C++ -*-

#ifndef __inline_myMeas_h__
#define __inline_myMeas_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"
#include "chroma_config.h"

namespace Chroma 
{ 
  // Environment in which the measurement lives (holds params and such)
  namespace InlineMyMeasIOGEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMyMeasIOGParams 
  {
    // Default constructor
    InlineMyMeasIOGParams();
    // Construct from XML
    InlineMyMeasIOGParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      int        cfg_serial;            /*!< The configuration serial number*/
      
      multi1d<std::string> hadrons;
      std::string l_prop;
      std::string file_name;

      // SftMomSrcPos_t is defined in lib/util/ft/sftmom.h
    } param;


    std::string xml_file; /*!< Alternate XML file pattern */
  }; // end of struct InlineMyMeasIOGParams


  //! Inline measurement of 3-pt functions writing building-blocks
  /*! \ingroup inlinehadron */
  class InlineMyMeasIOG : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineMyMeasIOG() {}
    // Constructor from param struct: copies param struct
    InlineMyMeasIOG(const InlineMyMeasIOGParams& p) : params(p) {}
    // Copy constructor
    InlineMyMeasIOG(const InlineMyMeasIOG& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineMyMeasIOGParams params;
  }; // end of class InlineMyMeasIOG






}; // end of namespace Chroma

#endif
