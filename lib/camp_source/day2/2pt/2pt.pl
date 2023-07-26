#!/usr/bin/perl
##

$conf=$ARGV[0];
$prefix=$ARGV[1];

$save_path="/dssg/home/acct-phyww/phyww/qazhang/training_camp/clqcd/C24P31/";
$quark_mass_l="-0.2770";
$quark_mass_s="-0.2400";

print <<"EOF";
<?xml version="1.0"?>
<chroma>
<Param> 
  <InlineMeasurements>

    <elem>
        <annotation>
            Read the L prop
        </annotation>
        <Name>QIO_READ_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>L_quark_propagator</object_id>
            <object_type>LatticePropagator</object_type>
        </NamedObject>
        <File>
            <file_name>${save_path}/Propagators/prop_${conf}_CoulombWall_t0-n1_P505050_m${quark_mass_l}_single</file_name>
            <parallel_io>true</parallel_io>
        </File>
    </elem>

    <elem>
        <annotation>
            Read the S prop
        </annotation>
        <Name>QIO_READ_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>S_quark_propagator</object_id>
            <object_type>LatticePropagator</object_type>
        </NamedObject>
        <File>
            <file_name>${save_path}/Propagators/prop_${conf}_CoulombWall_t0-n1_P505050_m${quark_mass_s}_single</file_name>
            <parallel_io>true</parallel_io>
        </File>
    </elem>


    <elem>
      <annotation>
         Compute the measurements you build
      </annotation>
      <Name>My_Measurements</Name>
      <Param>
        <cfg_serial>${conf}</cfg_serial>
        <hadrons>
          <elem>PION</elem>
          <elem>KAON</elem>
          <elem>PROTON</elem>
        </hadrons>
        <l_prop>L_quark_propagator</l_prop>
        <s_prop>S_quark_propagator</s_prop>
        <file_name>${prefix}/Data/2pt_${conf}.dat.iog</file_name>
       </Param>
    </elem>

  </InlineMeasurements>
    <nrow>24 24 24 72</nrow>
</Param>

  <RNG>
    <Seed>
      <elem>11</elem>
      <elem>11</elem>
      <elem>11</elem>
      <elem>0</elem>
    </Seed>
  </RNG>

  <Cfg>
    <cfg_type>SCIDAC</cfg_type>
    <cfg_file>${save_path}/Configurations/CoulombGaugeFixed/beta6.20_mu-0.2770_ms-0.2400_L24x72_cfg_${conf}_hyp0_gfixed3.scidac</cfg_file>
    <parallel_io>true</parallel_io>
  </Cfg>

</chroma>



EOF












