<?xml version="1.0" encoding="ISO-8859-1" ?>
<OFELI_File>
<info>
   <title></title>
   <date></date>
   <author></author>
</info>
<Project name="Step">
   <time_step value="0.005"/>
   <max_time value="5.0"/>
   <verbose value="1"/>
   <plot value="10"/>
   <mesh_file value="sstep.m"/>
   <parameter label="Reynolds" value="100"/>
   <parameter label="v_file" value="step.v"/>
   <parameter label="p_file" value="step.p"/>
</Project>
<Prescription>
   <BoundaryCondition code="2" dof="1">(y-1)*(3-y)</BoundaryCondition>
</Prescription>
</OFELI_File>
