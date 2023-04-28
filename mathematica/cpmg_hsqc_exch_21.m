(* SetDirectory["/Users/apple/Desktop"] *)
(* Remove["Global`@*"]   type this for each execution *)
(* << cpmg_hsqc_exch_21.m *)

(* bloch[1] shows peak intensity as a function of the offset *)
(* bloch[3] shows Rex as a function of the Vcpmg *)
(* bloch[2] simulation of 60-300-60 composite pulse against offset *)
(* bloch[4] simulation of 60-300-60 composite pulse against delay *)

hpi = 6;
(*
	5, ST-CW-CPMG {y y x -x}
	6, in-phase 2*cpmg(-Sy) - 180x - 2*cpmg(Sy) with 1H CW decoupling
	7, simple 2*cpmg(2SyIz) - Pelement - 2*cpmg(Sx), the conventional method
	8, AFTAC including the 60-300-60 pulses
*)

bloch[n_] :=
Module[
	{i, j, tmpi, tmps, tmpis, tmpisz, hamt},

	(* Hamiltonian for state A *)

	wxi = 2*Pi*v1i*Cos[phi];
	wyi = 2*Pi*v1i*Sin[phi];
	wzi = 2*Pi*offi;
	tmpi = KroneckerProduct[{{1, 0}, {0, 0}}, {{-r2i, -wzi, wyi}, {wzi, -r2i, -wxi}, {-wyi, wxi, -r1i}}];

	wxs = 2*Pi*v1s*Cos[phs];
	wys = 2*Pi*v1s*Sin[phs];
	wzs = 2*Pi*offs;
	tmps = KroneckerProduct[{{0, 0}, {0, 1}}, {{-r2s, -wzs, wys}, {wzs, -r2s, -wxs}, {-wys, wxs, -r1s}}];

	tmpis = tmpi + tmps;
	tmpi = KroneckerProduct[{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{-r2i, -wzi, wyi}, {wzi, -r2i, -wxi}, {-wyi, wxi, -r1i}}];
	tmps = KroneckerProduct[{{-r2s, -wzs, wys}, {wzs, -r2s, -wxs}, {-wys, wxs, -r1s}}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}];

	tmpisz = tmpi + tmps;
	tmpi = {{0,0,0,0,0,0,0,-Pi*jis,0}, {0,0,0,0,0,0,Pi*jis,0,0}, {0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,-Pi*jis,0,0,0}, {0,0,Pi*jis,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0}};
	tmps = Transpose[{{0,0,0,0,0,0,0,Pi*jis,0}, {0,0,0,0,0,0,-Pi*jis,0,0}, {0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,Pi*jis,0,0,0}, {0,0,-Pi*jis,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0}}];

	hamt = ArrayFlatten[{{tmpis, tmpi}, {tmps, tmpisz}}];

	(* Hamiltonian for state B *)

	bwzi = 2*Pi*(offi+dwib);
	bwzs = 2*Pi*(offs+dwsb);
	tmpi = KroneckerProduct[{{1, 0}, {0, 0}}, {{-r2i, -bwzi, wyi}, {bwzi, -r2i, -wxi}, {-wyi, wxi, -r1i}}];
	tmps = KroneckerProduct[{{0, 0}, {0, 1}}, {{-r2s, -bwzs, wys}, {bwzs, -r2s, -wxs}, {-wys, wxs, -r1s}}];

	tmpis = tmpi + tmps;
	tmpi = KroneckerProduct[{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{-r2i, -bwzi, wyi}, {bwzi, -r2i, -wxi}, {-wyi, wxi, -r1i}}];
	tmps = KroneckerProduct[{{-r2s, -bwzs, wys}, {bwzs, -r2s, -wxs}, {-wys, wxs, -r1s}}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}];

	tmpisz = tmpi + tmps;
	tmpi = {{0,0,0,0,0,0,0,-Pi*jis,0}, {0,0,0,0,0,0,Pi*jis,0,0}, {0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,-Pi*jis,0,0,0}, {0,0,Pi*jis,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0}};
	tmps = Transpose[{{0,0,0,0,0,0,0,Pi*jis,0}, {0,0,0,0,0,0,-Pi*jis,0,0}, {0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,Pi*jis,0,0,0}, {0,0,-Pi*jis,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0,0}}];

	hamtb = ArrayFlatten[{{tmpis, tmpi}, {tmps, tmpisz}}];

	(* Hamiltonian for state A *)

	hama = PadLeft[hamt, {16, 16}];

	(* auto- and cross-correlated cross-relaxation, ignored *)
	hama[[4, 7]] = 0;
	hama[[7, 4]] = 0;
	hama[[8, 12]] = 0;
	hama[[12, 8]] = 0;
	hama[[9, 11]] = 0;
	hama[[11, 9]] = 0;

	(* approximated relaxation *)
	hama[[8, 8]] = -r2s;	(* MQ *)
	hama[[9, 9]] = -r2s;	(* MQ *)
	hama[[10, 10]] = -r2s;	(* anti-phase, SxyIz *)
	hama[[11, 11]] = -r2s;	(* MQ *)
	hama[[12, 12]] = -r2s;	(* MQ *)
	hama[[13, 13]] = -r2s;	(* anti-phase, SxyIz *)
	hama[[14, 14]] = -r2i;	(* anti-phase, SzIxy *)
	hama[[15, 15]] = -r2i;	(* anti-phase, SzIxy *)
	hama[[16, 16]] = 0;		(* two-spin order *)

	(* Hamiltonian for state B *)

	(* auto- and cross-correlated cross-relaxation, ignored *)
	hamtb[[3, 6]] = 0;
	hamtb[[6, 3]] = 0;
	hamtb[[7, 11]] = 0;
	hamtb[[11, 7]] = 0;
	hamtb[[8, 10]] = 0;
	hamtb[[10, 8]] = 0;

	(* approximated relaxation *)
	hamtb[[7, 7]] = -r2s;	(* MQ *)
	hamtb[[8, 8]] = -r2s;	(* MQ *)
	hamtb[[9, 9]] = -r2s;	(* anti-phase, SxyIz *)
	hamtb[[10, 10]] = -r2s;	(* MQ *)
	hamtb[[11, 11]] = -r2s;	(* MQ *)
	hamtb[[12, 12]] = -r2s;	(* anti-phase, SxyIz *)
	hamtb[[13, 13]] = -r2i;	(* anti-phase, SzIxy *)
	hamtb[[14, 14]] = -r2i;	(* anti-phase, SzIxy *)
	hamtb[[15, 15]] = 0;	(* two-spin order *)

	ham = ArrayFlatten[{{hama, 0}, {0, hamtb}}];

	(* exchange involving auto- and cross-correlated cross-relaxation, ignored *)
	ham[[4, 1]] = 2 * r1i * kba/(kab + kba);
	ham[[7, 1]] = 2 * r1s * kba/(kab + kba);
	ham[[16, 1]] = 0;
	ham[[19, 1]] = 2 * r1i * kab/(kab + kba);
	ham[[22, 1]] = 2 * r1s * kab/(kab + kba);
	ham[[31, 1]] = 0;

	For
	[
		i = 2, i <= 16, i++,
		ham[[i, i]] -= kab
	];

	For
	[
		i = 17, i <= 31, i++,
		ham[[i, i]] -= kba
	];

	For
	[
		i = 2, i <= 16, i++,
		j = i + 15;
		ham[[i, j]] += kba;
		ham[[j, i]] += kab
	];

	vecini = {1/2, ix, iy, iz, sx, sy, sz, ixsx, iysx, izsx, ixsy, iysy, izsy, ixsz, iysz, izsz, bix, biy, biz, bsx, bsy, bsz, bixsx, biysx, bizsx, bixsy, biysy, bizsy, bixsz, biysz, bizsz};

	r2i = 0;
	r1i = 0;
	r2s = 3;
	r1s = 0;

	(* Do not input 0 to both kab and kba at the same time. *)
	kab = 200;	(* exchange rate from A to B *)
	kba = 200;	(* exchange rate from B to A *)

	jis = 140.0;	(* J(IS) coupling constant *)

	dwib = 0;	(* I's chemical shift difference of state B from that of state A *)
	dwsb = 50;	(* S's chemical shift difference of state B from that of state A *)

	phi = 0;
	v1i = 0;
	phs = 0;
	v1s = 0;
	offi = 0;
	offs = 3500;

	ix = 0;
	iy = 0;
	iz = 0;
	sx = 0;
	sz = 0;
	ixsx = 0;
	iysx = 0;
	izsx = 0;
	ixsy = 0;
	iysy = 0;
	ixsz = 0;
	iysz = 0;
	izsz = 0;

	bix = 0;
	biy = 0;
	biz = 0;
	bsx = 0;
	bsz = 0;
	bixsx = 0;
	biysx = 0;
	bizsx = 0;
	bixsy = 0;
	biysy = 0;
	bixsz = 0;
	biysz = 0;
	bizsz = 0;

	If
	[
		hpi == 7 || hpi == 8,
			izsy = -kba/(kab + kba);	(* the initial state *)
			bizsy = -kab/(kab + kba);
			sy = 0;
			bsy = 0,
		(* else hpi == 5, 6 *)
			izsy = 0;
			bizsy = 0;
			sy = -kba/(kab + kba);
			bsy = -kab/(kab + kba)
	];

	tex = 40*10^-3;		(* Total relaxation delay = 4 * each relax delay *)
	n90 = 20*10^-6;		(* 90 degree pulse width for spin S *)
	h90 = 10*10^-6;		(* 90 degree pulse width for spin I *)
	pmis = 0.9;			(* mis-calibration of the pulse width for spin S *)
	hpmis = 1.0;		(* mis-calibration of the pulse width for spin I *)
	dmis = 1.0;			(* mis-setting of 1/(2J) delay *)

	ofsrng = 4000.0;	(* The offset range of +- ofsrng, required for bloch[1] *)
	gmax = 40;			(* How many offset points are within +- ofsrng, required for bloch[1] *)
	ofs = offs;			(* The S-spin offset, required for bloch[3] *)
	lctm = 18;			(* How many Vcpmg points, required for bloch[3] *)
	minhcw = 15000;		(* Min power of 1H for CW decoupling *)

	Switch
	[
		n,
		1,	peakVSoffset[1],	(* peak intensity as a function of the offset *)
		3,	rexVSvcpmg[1],		(* Rex as a function of the Vcpmg *)
		2,	composite300[1],	(* simulation of 60-300-60 composite pulse against offset *)
		4,	comp300[1]			(* simulation of 60-300-60 composite pulse against delay *)
	]
];

(* Inphase 1H-CW pulse sequence *)
cwpulseInphase[de_, lct_, p1_, c1_] :=
Module[
	{k, pcw, hpwr, tphs1, cphs1},

	tphs1 = phs + p1 * Pi/2;
	cphs1 = phs + c1 * Pi/2;

	k = 1;
	If
	[
		lct == 0,
		pcw = tex/4.0,
		pcw = tex/lct/4.0		(* the length for one 13C echo = the length of a 1H 360deg CW pulse *)
	];
	While
	[
		(pcw/k) > 1/minhcw,
		k += 1
	];
	hpwr = k/pcw;

	(*
		The above code can adjust the 1H cw power such that Vcw = 2n * Vcpmg
		Activating the following line can fix the 1H cw power.
		# hpwr = minhcw;
	*)

	For
	[
		k = 1, k <= lct, k++,
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]]
	];

	For
	[
		k = 1, k <= lct, k++,
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]]
	];

	(* S: 60-300-60 *)

	vec = Block[{v1i = hpwr}, MatrixExp[ham*0.025/1000, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*300/90.0), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = Block[{v1i = hpwr}, MatrixExp[ham*0.025/1000, vec]];

	For
	[
		k = 1, k <= lct, k++,
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]]
	];

	For
	[
		k = 1, k <= lct, k++,
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]]
	]
];

(* 1H-CW yyx-x pulse sequence *)
cwpulseYyxx[de_, lct_] :=
Module[
	{k, pcw, hpwr},

	k = 1;
	If
	[
		lct == 0,
		pcw = tex/4.0,
		pcw = tex/lct/4.0		(* the length for one 13C echo = the length of a 1H 360deg CW pulse *)
	];
	While
	[
		(pcw/k) > 1/minhcw,
		k += 1
	];
	hpwr = k/pcw;

	(*
		The above code can adjust the 1H cw power such that Vcw = 2n * Vcpmg
		Activating the following line can fix the 1H cw power.
		# hpwr = minhcw;
	*)

	For
	[
		k = 1, k <= lct, k++,

		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = phs + 1 * Pi/2, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];

		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = phs + 1 * Pi/2, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];

		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = phs + 0 * Pi/2, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];

		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]];
		vec = Block[{v1s = 1/4/n90*pmis, phs = phs + 2 * Pi/2, v1i = hpwr}, MatrixExp[ham*n90*2, vec]];
		vec = Block[{v1i = hpwr}, MatrixExp[ham*de, vec]]
	]
];

(* The conventional RC pulse sequence *)
rcpulseConventional[de_, lct_, p1_, c1_, c2_, c3_, c4_] :=
Module[
	{k, tphs1, cphs1, cphs2, cphs3, cphs4},

	tphs1 = phs + p1 * Pi/2;
	cphs1 = phs + c1 * Pi/2;
	cphs2 = phs + c2 * Pi/2;
	cphs3 = phs + c3 * Pi/2;
	cphs4 = phs + c4 * Pi/2;

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs2}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	(* S: 60-300-60, I: 90-180-90 *)
	(* ******* n90*150/90.0-h90*2 must be > 0 ********)

	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = Pi/2}, MatrixExp[ham*h90*2, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs3}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs4}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	]
];

(* The new AFTAC pulse sequence *)
rcpulseAFTAC[de_, lct_, p1_, p2_, p3_, c1_, c2_, c3_, c4_] :=
Module[
	{k, tphs1, tphs2, tphs3, cphs1, cphs2, cphs3, cphs4},

	tphs1 = phs + p1 * Pi/2;
	tphs2 = phs + p2 * Pi/2;
	tphs3 = phs + p3 * Pi/2;
	cphs1 = phs + c1 * Pi/2;
	cphs2 = phs + c2 * Pi/2;
	cphs3 = phs + c3 * Pi/2;
	cphs4 = phs + c4 * Pi/2;

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs1}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	(* S: 60-300-60, I: 90-180-90 *)
	(* ******* n90*150/90.0-h90*2 must be > 0 ********)

	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = Pi/2}, MatrixExp[ham*h90*2, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs2}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	Switch
	[
		hpi,
		(* the central 180deg S pulse *)
			(* ?, vec = Block[{v1s = 1/4/n90*pmis, phs = tphs2}, MatrixExp[ham*n90*2, vec]], (* simple 180 deg *) *)
			8, vec = MatrixExp[dmis*ham*1.22/1000, vec];
			   vec = Block[{v1s = 1/4/n90*pmis, phs = tphs2+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
			   vec = Block[{v1s = 1/4/n90*pmis, phs = tphs2}, MatrixExp[ham*(n90*300/90.0), vec]];
			   vec = Block[{v1s = 1/4/n90*pmis, phs = tphs2+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
			   vec = MatrixExp[dmis*ham*1.22/1000, vec]
	];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs3}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	];

	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3, v1i = 1/4/h90*hpmis, phi = Pi/2}, MatrixExp[ham*h90*2, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
	vec = Block[{v1s = 1/4/n90*pmis, phs = tphs3+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
	vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*2.2), vec];

	For
	[
		k = 1, k <= lct, k++,
		vec = MatrixExp[ham*de, vec];
		vec = Block[{v1s = 1/4/n90*pmis, phs = cphs4}, MatrixExp[ham*n90*2, vec]];
		vec = MatrixExp[ham*de, vec]
	]
];

rexVSvcpmg[p_] :=
Module[
	{d17, lct},

	magx1 = Array[magxa1, {lctm, 2}];
	vec = Array[veca, {31, 1}];

	zg = OpenWrite["magmt.txt"];

	For
	[
		lct = 0, lct <= lctm, lct++,

		If
		[
			lct == 0,
			d17 = 0.0,
			d17 = tex/lct/8.0 - n90;
			magx1[[lct, 1]] = 2*lct/tex
		];

		vec = vecini;
		Switch
		[
			hpi,
			5, Block[{offs = ofs}, cwpulseYyxx[d17, lct]]; (* 1H CW yyxx *)
				rec1 = vec[[6]],
			6, Block[{offs = ofs}, cwpulseInphase[d17, lct, 0, 1]]; (* 1H CW inphase 2*cpmg(-Sy) - 180x - 2*cpmg(Sy) *)
				rec1 = vec[[6]],
			7, Block[{offs = ofs}, rcpulseConventional[d17, lct, 0, 1, 1, 0, 0]]; (* 2*cpmg(2SyIz) - P - 2*cpmg(Sx) *)
				rec1 = vec[[5]],
			8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 0, 2, 1, 0, 2, 1]]; (* 60-300-60 pulse *)
				rec1 = vec[[13]]
		];

		vec = vecini;
		Switch
		[
			hpi,
			5, rec2 = 0,
			6, Block[{offs = ofs}, cwpulseInphase[d17, lct, 2, 1]]; (* 1H CW inphase 2*cpmg(-Sy) - 180-x - 2*cpmg(Sy) *)
				rec2 = vec[[6]],
			7, Block[{offs = ofs}, rcpulseConventional[d17, lct, 2, 1, 1, 0, 0]]; (* 2*cpmg(2SyIz) - P - 2*cpmg(Sx) *)
				rec2 = vec[[5]],
			8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 2, 2, 1, 0, 2, 1]]; (* 60-300-60 pulse *)
				rec2 = vec[[13]]
		];

		rec3 = 0;
		rec4 = 0;

		Switch
		[
			hpi,
			5, ns = 1,
			6, ns = 2,
			7, ns = 2,
			8, ns = 2
		];
		yrec = N[(rec1 + rec2 + rec3 + rec4)/ns];

		If
		[
			lct == 0,
			iyrec = yrec,
			magx1[[lct, 2]] = -Log[N[yrec/iyrec]] / tex;
			WriteString[zg, magx1[[lct, 1]], "	", magx1[[lct, 2]], "	", magx1[[1, 2]]*0.02,"\n"]
		]
	];
	Close[zg];

	g1 = ListPlot[{magx1}, Joined->True, Mesh->Full, PlotLegends->{"standard"}, Frame -> True, PlotStyle->Cyan, PlotRange -> {{0, 1050}, {-1, 150}}]
];

peakVSoffset[p_] :=
Module[
	{h, g, d17, lct, ofs, ns},

	magx = Array[magxa, {gmax, 2}];
	magx0 = Array[magxa0, {gmax, 2}];
	magx1 = Array[magxa1, {gmax, 2}];
	magx2 = Array[magxa2, {gmax, 2}];
	magx3 = Array[magxa3, {gmax, 2}];
	vec = Array[veca, {31, 1}];

	For
	[
		h = 0, h <= 3, h++,

		Switch
		[
			h,
			0,	lct = 0,
			1,	lct = 1,
			2,	lct = 2,
			3,	lct = 40
		];

		If
		[
			lct == 0, d17 = 0, d17 = tex/lct/8.0 - n90
		]

		For
		[
			g = 1, g <= gmax, g++,

			ofs = (g-gmax/2.0)*2*ofsrng/gmax;
			rec1 = 0;
			rec2 = 0;
			rec3 = 0;
			rec4 = 0;

			vec = vecini;
			Switch
			[
				hpi,
				5, Block[{offs = ofs}, cwpulseYyxx[d17, lct]]; (* 1H CW yyxx *)
					rec1 = vec[[6]],
				6, Block[{offs = ofs}, cwpulseInphase[d17, lct, 0, 1]]; (* 1H CW inphase 2*cpmg(-Sy) - 180x - 2*cpmg(Sy) *)
					rec1 = vec[[6]],
				7, Block[{offs = ofs}, rcpulseConventional[d17, lct, 0, 1, 1, 0, 0]]; (* 2*cpmg(2SyIz) - P - 2*cpmg(Sx) *)
					rec1 = vec[[5]],
					(*
					rec1 = (vec[[5]] - vec[[10]])*0.5,	(* TROSY component *)
					rec1 = (vec[[5]] + vec[[10]])*0.5,	(* anti-TROSY component *)
					*)

				8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 0, 2, 1, 0, 2, 1]]; (* 60-300-60 pulse Best *)
				(*
				8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 0, 2, 3, 0, 2, 3]]; (* 60-300-60 Varian *)
				*)
					rec1 = vec[[13]]	(* state A *)
				(*
					rec1 = vec[[28]]	(* +15, state B *)
				*)
			];

			vec = vecini;
			Switch
			[
				hpi,
				5, rec2 = 0,
				6, Block[{offs = ofs}, cwpulseInphase[d17, lct, 2, 1]]; (* 1H CW inphase 2*cpmg(-Sy) - 180-x - 2*cpmg(Sy) *)
					rec2 = vec[[6]],
				7, Block[{offs = ofs}, rcpulseConventional[d17, lct, 2, 1, 1, 0, 0]]; (* 2*cpmg(2SyIz) - P - 2*cpmg(Sx) *)
					rec2 = vec[[5]],
					(*
					rec2 = (vec[[5]] - vec[[10]])*0.5,	(* TROSY component *)
					rec2 = (vec[[5]] + vec[[10]])*0.5,	(* anti-TROSY component *)
					*)

				8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 2, 2, 1, 0, 2, 1]]; (* 60-300-60 pulse Best *)
				(*
				8, Block[{offs = ofs}, rcpulseAFTAC[d17, lct, 0, 2, 2, 3, 0, 2, 3]]; (* 60-300-60 pulse Varian *)
				*)
					rec2 = vec[[13]]	(* state A *)
				(*
					rec2 = vec[[28]]	(* +15, state B *)
				*)
			];

			Switch
			[
				hpi,
				5, ns = 1,
				6, ns = 2,
				7, ns = 2,
				8, ns = 2
			];

			(*
			magx[[g, 2]] = N[(rec1 + rec2 + rec3 + rec4)/ns/(-vecini[[28]])];	(* state B *)
			magx[[g, 2]] = N[(rec1 + rec2 + rec3 + rec4)/ns/(-vecini[[13]])];	(* state A *)
			*)
			magx[[g, 2]] = N[(rec1 + rec2 + rec3 + rec4)/ns];
			magx[[g, 1]] = ofs
		];

		Switch
		[
			h,
			0,	magx0 = magx,
			1,	magx1[[All, 2]] = magx[[All, 2]] / magx0[[All, 2]] / 1.0; magx1[[All, 1]] = magx[[All, 1]],
			2,	magx2[[All, 2]] = magx[[All, 2]] / magx0[[All, 2]] / 1.0; magx2[[All, 1]] = magx[[All, 1]],
			3,	magx3[[All, 2]] = magx[[All, 2]] / magx0[[All, 2]] / 1.0; magx3[[All, 1]] = magx[[All, 1]]
		]
	];

	gg = ListPlot[{magx1, magx2, magx3},
	Joined -> {True},
	PlotLegends -> {"n=1", "n=2", "n=40"},
	(* PlotStyle -> {Directive[Black], Directive[Black, Dashed], Directive[Gray, Thick]}, *)
	PlotStyle -> {Directive[Magenta, Thick], Directive[Blue, Dashed, Thick], Directive[Green, Dotted, Thick]},
	Axes -> False,
	Frame -> True,
	FrameStyle -> AbsoluteThickness[2.0],
	PlotRange -> {{-ofsrng, ofsrng}, {0, 1.1}},
	GridLines -> None,
	FrameTicksStyle -> Directive[FontFamily -> "Times New Roman", FontSize -> 17]
	]

	(*
	tofile["temp1.txt", magx1];
	tofile["temp2.txt", magx2];
	tofile["temp3.txt", magx3];
	Export["temp1.pdf", gg]
	*)
];

composite300[p_] :=
Module[
	{g, ofs},

	magx1 = Array[magxa1, {gmax, 2}];
	magx2 = Array[magxa2, {gmax, 2}];
	vec = Array[veca, {31, 1}];

	For
	[
		g = 1, g <= gmax, g++,

		ofs = (g-gmax/2.0)*2*ofsrng/gmax;
		vec = vecini;
		Block[{offs = ofs}, rcpulse300[1, 1]]; (* 60-300-60 pulse *)

		magx1[[g, 2]] = vec[[5]];	(* state A *)
		magx1[[g, 1]] = ofs;
		magx2[[g, 2]] = vec[[13]];
		magx2[[g, 1]] = ofs
	];

	(*
	vecini = {1/2, ix, iy, iz, sx, sy, sz, ixsx, iysx, izsx, ixsy, iysy, izsy, ixsz, iysz, izsz, bix, biy, biz, bsx, bsy, bsz, bixsx, biysx, bizsx, bixsy, biysy, bizsy, bixsz, biysz, bizsz};
	*)

	gg = ListPlot[{magx1, magx2},
	Joined -> {True},
	PlotLegends -> {"Sx", "IzSy"},
	(* PlotStyle -> {Directive[Black], Directive[Black, Dashed], Directive[Gray, Thick]}, *)
	PlotStyle -> {Directive[Magenta, Thick], Directive[Blue, Thick]},
	Axes -> False,
	Frame -> True,
	FrameStyle -> AbsoluteThickness[2.0],
	PlotRange -> {{-ofsrng, ofsrng}, {-0.5, 0.2}},
	GridLines -> None,
	FrameTicksStyle -> Directive[FontFamily -> "Times New Roman", FontSize -> 17]
	]

	(*
	tofile["temp1.txt", magx1];
	tofile["temp2.txt", magx2];
	Export["temp1.pdf", gg]
	*)
];

rcpulse300[p1_, pp_, dd_] :=
Module[
	{tphs1},

	tphs1 = phs + p1 * Pi/2;

	Switch
	[
		pp,
		1,
			(* 60-300-60 *)
			(* dd = 2.2 is the best (1.9-2.6), depending on the offset. *)

			vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*dd), vec];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = Pi/2}, MatrixExp[ham*h90*2, vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1, v1i = 1/4/h90*hpmis, phi = 0}, MatrixExp[ham*h90, vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1}, MatrixExp[ham*(n90*150/90.0-h90*2), vec]];
			vec = Block[{v1s = 1/4/n90*pmis, phs = tphs1+Pi}, MatrixExp[ham*n90*60/90.0, vec]];
			vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]-n90*2.0/Pi*dd), vec],

		2,
			(* normal 180 deg *)
			(* ****** The pulse length of S should be just twice that of I *****)

			hh90 = 10*10^-6;
			nn90 = 2 * hh90;

			vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]), vec];
			vec = Block[{v1s = 1/4/nn90*pmis, phs = tphs1, v1i = 1/4/hh90*hpmis, phi = 0}, MatrixExp[ham*hh90, vec]];
			vec = Block[{v1s = 1/4/nn90*pmis, phs = tphs1, v1i = 1/4/hh90*hpmis, phi = Pi/2}, MatrixExp[ham*hh90*2, vec]];
			vec = Block[{v1s = 1/4/nn90*pmis, phs = tphs1, v1i = 1/4/hh90*hpmis, phi = 0}, MatrixExp[ham*hh90, vec]];
			vec = MatrixExp[dmis*ham*(1.0/4/Abs[jis]), vec]
	]
];

(* Customize the 60-300-60 pulse *)
comp300[p_] :=
Module[
	{k, ofs},

	magx1 = Array[magxa1, {100, 2}];
	magx2 = Array[magxa2, {100, 2}];
	(* vec = Array[veca, {31, 1}]; *)

	ofs = 0;
	dd = 0;
	For
	[
		k = 1, k <= 100, k++,

		dd = k*0.1;
		vec = vecini;
		Block[{offs = ofs}, rcpulse300[0, 1, dd]]; (* 60-300-60 pulse *)

		magx1[[k, 2]] = vec[[5]];	(* state A *)
		magx1[[k, 1]] = dd;
		magx2[[k, 2]] = vec[[13]];
		magx2[[k, 1]] = dd
	];

	(*
	vecini = {1/2, ix, iy, iz, sx, sy, sz, ixsx, iysx, izsx, ixsy, iysy, izsy, ixsz, iysz, izsz, bix, biy, biz, bsx, bsy, bsz, bixsx, biysx, bizsx, bixsy, biysy, bizsy, bixsz, biysz, bizsz};
	*)

	gg = ListPlot[{magx1, magx2},
	Joined -> {True},
	PlotLegends -> {"Sx", "IzSy"},
	PlotStyle -> {Directive[Magenta, Thick], Directive[Blue, Thick]},
	Axes -> False,
	Frame -> True,
	FrameStyle -> AbsoluteThickness[2.0],
	PlotRange -> {{0, 10.1}, {-1.1, 1.1}},
	FrameTicksStyle -> Directive[FontFamily -> "Times New Roman", FontSize -> 17],
	GridLines -> Automatic
	]

	(*
	tofile["temp1.txt", magx1];
	tofile["temp2.txt", magx2];
	Export["temp1.pdf", gg]
	*)
]
