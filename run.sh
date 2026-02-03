#!/bin/bash
cd build
CONFID=0
for BETA in 2 4 6; do
	if [ "$BETA" -eq 2 ]; then
        Mc="-018"
        Massc="-01884"
        mc=-0.18840579710144945
		M1="-1868"
        Mass1="-01868"
        m1=-0.1868
		M2="-0968"
        Mass2="-00968"
        m2=-0.0968
	elif [ "$BETA" -eq 4 ]; then
        Mc="-01023"
        Massc="-01023"
        mc=-0.1023
		M1="-0933"
        Mass1="-00933"
        m1=-0.0933
		M2="-0033"
        Mass2="-00033"
        m2=-0.0033
	elif [ "$BETA" -eq 6 ]; then
        Mc="-0709"
        Massc="-00709"
        mc=-0.0709
		M1="-0619"
        Mass1="-00619"
        m1=-0.0619
		M2="-0281"
        Mass2="-00281"
        m2=-0.0281
	fi
	for N in 128; do
        echo "================== BETA = ${BETA}, N = ${N}, M0 = ${mc} =================="
        echo ""
	    printf "${BETA}\n${mc}\n${CONFID}\n../../SchwingerModel/fermions/SchwingerModel/confs/b${BETA}_${N}x${N}/m${Mc}/2D_U1_Ns${N}_Nt${N}_b${BETA}0000_m${Massc}_${CONFID}.ctxt\n../../SchwingerModel/fermions/SchwingerModel/confs/rhs/rhs_conf${CONFID}_${N}x${N}_b${BETA}0000_m${Massc}.rhs" >> inputs
        ./AMG_${N}x${N}_l4 < inputs
        rm inputs
        echo ""
        echo ""

        echo "================== BETA = ${BETA}, N = ${N}, M0 = ${m1} =================="
        echo ""
	    printf "${BETA}\n${m1}\n${CONFID}\n../../SchwingerModel/fermions/SchwingerModel/confs/b${BETA}_${N}x${N}/m${M1}/2D_U1_Ns${N}_Nt${N}_b${BETA}0000_m${Mass1}_${CONFID}.ctxt\n../../SchwingerModel/fermions/SchwingerModel/confs/rhs/rhs_conf${CONFID}_${N}x${N}_b${BETA}0000_m${Mass1}.rhs" >> inputs
        ./AMG_${N}x${N}_l4 < inputs
        rm inputs
        echo ""
        echo ""
        
        echo "================== BETA = ${BETA}, N = ${N}, M0 = ${m2} =================="
        echo ""
	    printf "${BETA}\n${m2}\n${CONFID}\n../../SchwingerModel/fermions/SchwingerModel/confs/b${BETA}_${N}x${N}/m${M2}/2D_U1_Ns${N}_Nt${N}_b${BETA}0000_m${Mass2}_${CONFID}.ctxt\n../../SchwingerModel/fermions/SchwingerModel/confs/rhs/rhs_conf${CONFID}_${N}x${N}_b${BETA}0000_m${Mass2}.rhs" >> inputs
        ./AMG_${N}x${N}_l4 < inputs
        rm inputs
        echo ""
        echo ""
    done
done
