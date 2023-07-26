// $Id: sh_source_const.cc,v 3.14 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "grid_source.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/quark_displacement_factory.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/no_quark_displacement.h"

#include "meas/sources/zN_src.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const std::string& path, GridQuarkSourceConstEnv::Params& param)
  {
    GridQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const GridQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Hooks to register the class
  namespace GridQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }

      //! Callback function
      QuarkSourceConstruction<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
									  const std::string& path)
      {
	return new SourceConst<LatticeStaggeredPropagator>(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name("MOM_GRID_SOURCE");
    }

    //! Return the name
    std::string getName() {return std::string("SHELL_SOURCE");}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LinkSmearingEnv::registerAll();
	success &= QuarkSmearingEnv::registerAll();
	success &= QuarkDisplacementEnv::registerAll();
	success &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSourceConstructionFactory::Instance().registerObject(name, createStagProp);

	registered = true;
      }
      return success;
    }


    //! Read parameters
    Params::Params()
    {
      j_decay = -1;
      t_srce.resize(Nd);
      t_srce = 0;
      quark_smear_lastP = true;
    }

    //! Read parameters
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      quark_smear_lastP = true;

      switch (version) 
      {
      case 1:
      {
	quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
	quark_displacement.id = SimpleQuarkDisplacementEnv::getName();
	int disp_length = 0;
	int disp_dir = 0;

	XMLBufferWriter xml_tmp;
	push(xml_tmp, "Displacement");
	write(xml_tmp, "DisplacementType", quark_displacement.id);

	if (paramtop.count("disp_length") != 0)
	  read(paramtop, "disp_length", disp_length);

	if (paramtop.count("disp_dir") != 0)
	  read(paramtop, "disp_dir", disp_dir);

	write(xml_tmp, "disp_length", disp_length);
	write(xml_tmp, "disp_dir",  disp_dir);

	pop(xml_tmp);  // Displacement

	quark_displacement.xml = xml_tmp.printCurrentContext();
      }
      break;

      case 2:
      {
	if (paramtop.count("Displacement") != 0)
	  quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
	else
	  quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
      }
      break;

      case 3:
      {
	read(paramtop, "quark_smear_lastP", quark_smear_lastP);

	if (paramtop.count("Displacement") != 0)
	  quark_displacement = readXMLGroup(paramtop, "Displacement", "DisplacementType");
	else
	  quark_displacement = QuarkDisplacementEnv::nullXMLGroup();
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }

      quark_smearing = readXMLGroup(paramtop, "SmearingParam", "wvf_kind");

      if (paramtop.count("LinkSmearing") != 0)
	link_smearing = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      else
	link_smearing = LinkSmearingEnv::nullXMLGroup();

      read(paramtop, "t_srce", t_srce);
      read(paramtop, "grid", t_grid);
      read(paramtop, "j_decay",  j_decay);
      read(paramtop, "ini_mom",  ini_mom);

      Z3_phase=false;
      if (paramtop.count("Z3_phase") != 0)
      {
         read(paramtop, "Z3_phase",  Z3_phase);
      }
      
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 3;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", GridQuarkSourceConstEnv::name);
      xml << quark_smearing.xml;
      xml << quark_displacement.xml;
      xml << link_smearing.xml;

      write(xml, "t_srce",  t_srce);
      write(xml, "grid",  t_grid);
      write(xml, "j_decay",  j_decay);
      write(xml, "ini_mom", ini_mom);
      write(xml, "quark_smear_lastP",  quark_smear_lastP);

      pop(xml);
    }
    
    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Shell source" << std::endl;

      LatticePropagator quark_source;

      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	{
	  std::istringstream  xml_l(params.link_smearing.xml);
	  XMLReader  linktop(xml_l);
	  QDPIO::cout << "Link smearing type = " << params.link_smearing.id << std::endl;
	
	  Handle< LinkSmearing >
	    linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.link_smearing.id,
									 linktop,
									 params.link_smearing.path));
	  (*linkSmearing)(u_smr);
	}

	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
        QDPIO::cout << "Quark smearing type = " << params.quark_smearing.id << std::endl;
	
	Handle< QuarkSmearing<LatticePropagator> >
	  quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing.id,
									smeartop,
									params.quark_smearing.path));

	//
	// Create the quark displacement object
	//
	std::istringstream  xml_d(params.quark_displacement.xml);
	XMLReader  displacetop(xml_d);
        QDPIO::cout << "Displacement type = " << params.quark_displacement.id << std::endl;
	
	Handle< QuarkDisplacement<LatticePropagator> >
	  quarkDisplacement(ThePropDisplacementFactory::Instance().createObject(params.quark_displacement.id,
										displacetop,
										params.quark_displacement.path));


/// set the phase.

    LatticeReal p_dot_x;
    p_dot_x = 0.;
    for (int mu = 0, j=0; mu < Nd; ++mu) {
       const Real twopi = 6.283185307179586476925286;
       if (mu ==params.j_decay) continue;
        
       p_dot_x += LatticeReal(Layout::latticeCoordinate(mu) ) * twopi *
            Real(params.ini_mom[j]) / Layout::lattSize()[mu];
        
       ++j;
    }
    LatticeComplex  phase=cmplx(cos(p_dot_x),sin(p_dot_x));

/// set the phase.
	//
	// Create quark source
	//
	LatticeComplex c;
	LatticeReal rnd1, theta;
	random(rnd1);
	Real twopiN = Chroma::twopi / 3; //Z_3 grid for bayron;
	if(params.Z3_phase==true)
     	   theta = twopiN * floor(3*rnd1);
        else
           theta = zero;
	c = cmplx(cos(theta),sin(theta));

	multi1d<int> offset(4);
	for(int idr=0;idr<4;idr++) offset[idr]=(params.t_grid[idr]>0)?params.t_grid[idr]:(Layout::lattSize()[idr]);
	QDPIO::cout << "Offsets: " << offset[0] << " , " << offset[1] << " , " << offset[2] 
	      << " , " << offset[3] <<   std::endl;
	LatticeComplex c1=where(
	  ((((Layout::latticeCoordinate(0))-params.t_srce[0])%offset[0]==0)&&
	   (((Layout::latticeCoordinate(1))-params.t_srce[1])%offset[1]==0)&&
	   (((Layout::latticeCoordinate(2))-params.t_srce[2])%offset[2]==0)&&
	   (((Layout::latticeCoordinate(3))-params.t_srce[3])%offset[3]==0)),
	  c, LatticeComplex(zero));
       c1=c1*phase;	  
       for(int color_source = 0; color_source < Nc; ++color_source)
       for(int spin_source = 0; spin_source < Ns; ++spin_source)
       {
           LatticeColorVector colorvec = zero;
           LatticeFermion chi = zero; 
           pokeSpin(chi,pokeColor(colorvec,c1,color_source),spin_source);
           FermToProp(chi, quark_source, color_source, spin_source);
       }

	// Smear and displace
	if (params.quark_smear_lastP)
	{
	  // Smear the colour source
	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);

	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);
	}
	else
	{
	  // do the smearing
	  (*quarkSmearing)(quark_source, u_smr);

	  // Smear the colour source
	  // smear the point source first, then displace
	  // displacement has to be taken along negative direction.
	  (*quarkDisplacement)(quark_source, u_smr, MINUS);
	}

      }
      catch(const std::string& e) 
      {
        QDPIO::cerr << name << ": Caught Exception smearing: " << e << std::endl;
        QDP_abort(1);
      }

      LatticeComplex corr=trace(adj(quark_source)*quark_source);
      double value=sum(corr).elem().elem().elem().real();
      QDPIO::cout << "Norm2 of the source is " << value << std::endl;

      return quark_source;
    }

    //! Construct the source
    template<>
    LatticeStaggeredPropagator
    SourceConst<LatticeStaggeredPropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cerr << name << ": Caught Exception grid source: Staggered Fermion is not supported  " << std::endl;
      QDP_abort(1);
    }
  }
}
