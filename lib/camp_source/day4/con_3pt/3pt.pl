#!/usr/bin/perl
##

$conf=$ARGV[0];
$prefix=$ARGV[1];


$save_path="../../class4_prop/C24P31";
$l_mass="-0.2400";
$seqp_path="../Save_prop_data/Save_SeqProp";



print <<"EOF";
<?xml version="1.0"?>
<chroma>
<Param> 
  <InlineMeasurements>
EOF

@tseq_list = (3,4,5,6,7);
foreach $tseq (@tseq_list){
print <<"EOF";
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
            <file_name>${save_path}/Propagators/prop_${conf}_CoulombWall_t0-n1_P505050_m${l_mass}_single</file_name>
            <parallel_io>true</parallel_io>
        </File>
    </elem>

    <elem>
        <annotation>
            Read the Seq prop
        </annotation>
        <Name>QIO_READ_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>Seq_propagator</object_id>
            <object_type>LatticePropagator</object_type>
        </NamedObject>
        <File>
            <file_name>${seqp_path}/seqprop_P0_m${l_mass}_tseq${tseq}.${conf}</file_name>
            <parallel_io>true</parallel_io>
        </File>
    </elem>

    <elem>
      <annotation>
         Compute the measurements you build
      </annotation>
      <Name>Measure_gAgV</Name>
      <Param>
        <cfg_serial>${conf}</cfg_serial>
        <gAgV_curr>
          <elem>gV</elem>
          <elem>gA</elem>
        </gAgV_curr>
        <l_prop>L_quark_propagator</l_prop>
        <seq_prop>Seq_propagator</seq_prop>
        <file_name>${prefix}/Data/3pt_${conf}_tseq${tseq}.dat.iog</file_name>
       </Param>
    </elem>

    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>L_quark_propagator</object_id>
        </NamedObject>
    </elem>
    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>Seq_propagator</object_id>
        </NamedObject>
    </elem>

EOF
}

print <<"EOF";
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












