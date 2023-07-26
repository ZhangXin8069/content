#!/usr/bin/perl
##

$conf=$ARGV[0];
$prefix=$ARGV[1];

$save_path="../../class4_prop/C24P31";
$l_mass="-0.2400";

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
	    <file_name>${save_path}/Propagators/prop_${conf}_CoulombWall_t0-n1_P505050_m${l_mass}_single</file_name>
            <parallel_io>true</parallel_io>
        </File>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>1</version>
          <SinkType>POINT_SINK</SinkType>
          <j_decay>3</j_decay>
        </Sink>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>L_quark_propagator</prop_id>
        <smeared_prop_id>smeared_L_quark_propagator</smeared_prop_id>
      </NamedObject>
    </elem>
    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
            <object_id>L_quark_propagator</object_id>
        </NamedObject>
    </elem>

EOF

@tseq_list = (3,4,5,6,7);
foreach $tseq (@tseq_list){

print <<"EOF";
    <elem>
      <Name>SEQSOURCE_FAST</Name>
      <SmearedProps>true</SmearedProps>
      <multi_tSinks>${tseq}</multi_tSinks>
      <Frequency>1</Frequency>
      <Param>
        <version>2</version>
        <SeqSource>
          <version>1</version>
          <SeqSourceType>WEAK_CURRENT_gAgV_CUR</SeqSourceType>
          <j_decay>3</j_decay>
          <t_sink>0</t_sink>
          <sink_mom>0 0 0</sink_mom>
        </SeqSource>
      </Param>
      <PropSink>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>POINT_SINK</SinkType>
          <j_decay>3</j_decay>
          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>
        </Sink>
      </PropSink>
      <NamedObject>
        <prop_ids>
          <elem>smeared_L_quark_propagator</elem>
          <elem>smeared_L_quark_propagator</elem>
        </prop_ids>
        <seqsource_id>
          <elem>seqsrc_tseq_${tseq}</elem>
        </seqsource_id>
        <gauge_id>default_gauge_field</gauge_id>
      </NamedObject>
   </elem>

  <elem>
    <Name>PROPAGATOR</Name>
    <Frequency>1</Frequency>
    <Param>
      <version>10</version>
      <quarkSpinType>FULL</quarkSpinType>
      <obsvP>true</obsvP>
      <numRetries>1</numRetries>
        <FermionAction>
          <FermAct>UNPRECONDITIONED_CLOVER</FermAct>
          <Mass>${l_mass}</Mass>
          <clovCoeff>1.160920226</clovCoeff>
          <FermState>
              <Name>STOUT_FERM_STATE</Name>
              <rho>0.125</rho>
              <n_smear>1</n_smear>
              <orthog_dir>-1</orthog_dir>
              <FermionBC>
                <FermBC>SIMPLE_FERMBC</FermBC>
                <boundary>1 1 1 -1</boundary>
              </FermionBC>
            </FermState>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>1.0e-5</RsdCG>
          <MaxCG>1000</MaxCG>
        </InvertParam>
    </Param>
    <NamedObject>
      <gauge_id>default_gauge_field</gauge_id>
      <source_id>seqsrc_tseq_${tseq}</source_id>
      <prop_id>prop_P0_tseq_${tseq}</prop_id>
    </NamedObject>
  </elem>


    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
        <object_id>seqsrc_tseq_${tseq}</object_id>
        </NamedObject>
    </elem>

  <elem>
    <annotation>
      Write the named object
    </annotation>
    <Name>QIO_WRITE_NAMED_OBJECT</Name>
    <Frequency>1</Frequency>
    <NamedObject>
      <object_id>prop_P0_tseq_${tseq}</object_id>
      <object_type>LatticePropagator</object_type>
    </NamedObject>
    <File>
      <file_name>$prefix/seqprop_P0_m${l_mass}_tseq${tseq}.${conf}</file_name>
      <file_volfmt>SINGLEFILE</file_volfmt>
      <parallel_io>true</parallel_io>
    </File>
  </elem>


    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
        <object_id>prop_P0_tseq_${tseq}</object_id>
        </NamedObject>
    </elem>

EOF

}

print <<"EOF";
    <elem>
        <Name>ERASE_NAMED_OBJECT</Name>
        <Frequency>1</Frequency>
        <NamedObject>
        <object_id>smeared_L_quark_propagator</object_id>
        </NamedObject>
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












