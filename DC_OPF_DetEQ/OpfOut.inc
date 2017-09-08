FILE OPF_RES/ OPF_results.txt /;
PUT OPF_RES;

PUT 'Build Line l';
PUT /
LOOP(l,
	LOOP(i,
   		LOOP(j $ (y.l(l,i,j)=1),
			PUT 'Line ' l.TL:0;
			PUT y.l(l,i,j);
			PUT /;
		);
	);
);

PARAMETERS
xi(i,subs)            Real Power generated in node i
;

xi(i,subs) = SUM[g $ ig(i,g), x.l(g,subs)];

PUT 'Generated power [MW]':30;
PUT 'Generation':12;
PUT /;
LOOP(subs,
	LOOP(g,
   		PUT g.TL:30;
   		PUT x.l(g,subs):12;
   		PUT /;
	);
);

PUT /;
PUT /;

PUT 'Power in each node [GW]':30;
PUT 'Real Generation ':16;
PUT 'Load ':16;
PUT 'Angle ':16;
PUT 'Marg. cost ':16;
PUT 'Non-served':16;
PUT /;
LOOP(subs,
	PUT subs.TL:0;
	PUT /;
	LOOP(i,
		PUT i.TL:30;
		PUT xi(i,subs):12:4;
		PUT pLoad(i,subs):16:4;
		PUT theta.l(i,subs):16:4;
		PUT eDemand.m(i,subs):16:4;
		PUT pns.l(i,subs):16:4;
		PUT /;
	);
);

PUT /;
PUT /;

PUT 'Real Power flows between nodes [GW]':30;
PUT /;
LOOP(l,
	LOOP(i,
   		LOOP(j $ ii(l,i,j),
      		PUT 'Trans line from ' i.TL:0 ' to ' PUT j.TL:0 ' on ' PUT l.TL:0
			LOOP(subs,
		      		PUT p.l(l,i,j,subs):12:4;
			);
      		PUT /;
   		);
	);
);


$ontext
PUT /;
PUT /;
PUT 'Power losses [GW]               aprox     real':60;
PUT /;
LOOP(i,
   LOOP(j $ ii(i,j),
      PUT 'Trans line from ' i.TL:0 ' to ' PUT j.TL:0
      PUT l.l(i,j):14:5;
      PUT (2 * [dtl(i,j,'r')/ (dtl(i,j,'r')**2 + dtl(i,j,'x')**2)]* [1-COS(theta.l(i) - theta.l(j))]):12:5
      PUT /;
   );
);
$offtext
