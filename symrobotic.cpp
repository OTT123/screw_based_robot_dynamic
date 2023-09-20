#include "symrobotic.h"

Eigen::Matrix<double, DOF* DOF, 1> M_gluon(const double* parms, const double* q)
{
	Eigen::Matrix<double, DOF* DOF, 1> M;
	double x0 = cos(q[1]);
	double x1 = sin(q[1]);
	double x2 = -x1;
	double x3 = cos(q[2]);
	double x4 = -x0;
	double x5 = sin(q[2]);
	double x6 = x2 * x3 + x4 * x5;
	double x7 = -x6;
	double x8 = cos(q[3]);
	double x9 = x0 * x3 + x2 * x5;
	double x10 = sin(q[3]);
	double x11 = x10 * x9 + x7 * x8;
	double x12 = x10 * x6 + x8 * x9;
	double x13 = -0.08 * x1;
	double x14 = -x13;
	double x15 = -0.08 * x0;
	double x16 = x14 * x5 + x15 * x3;
	double x17 = -x9;
	double x18 = x16 + 0.005 * x17;
	double x19 = x13 * x3 + x15 * x5;
	double x20 = x19 + 0.005 * x6;
	double x21 = x10 * x18 + x20 * x8;
	double x22 = -x21;
	double x23 = -0.175 * x2;
	double x24 = x23 + 0.17 * x7;
	double x25 = cos(q[4]);
	double x26 = x11 * x25;
	double x27 = -x26;
	double x28 = sin(q[5]);
	double x29 = cos(q[5]);
	double x30 = x12 * x29 + x27 * x28;
	double x31 = x12 * x28 + x26 * x29;
	double x32 = -x11;
	double x33 = sin(q[4]);
	double x34 = x32 * x33;
	double x35 = -x34;
	double x36 = -0.08 * x11 - x24;
	double x37 = -x8;
	double x38 = x10 * x20 + x18 * x37;
	double x39 = x25 * x38 + x33 * x36;
	double x40 = x21 * x28 + x29 * x39;
	double x41 = -parms[56];
	double x42 = -x38;
	double x43 = x25 * x36 + x33 * x42;
	double x44 = -x43;
	double x45 = parms[51] * x31 + parms[53] * x30 + parms[54] * x35 + parms[58] * x40 + x41 * x44;
	double x46 = -x28;
	double x47 = -parms[58];
	double x48 = -x39;
	double x49 = x21 * x29 + x28 * x48;
	double x50 = parms[50] * x31 + parms[51] * x30 + parms[52] * x35 + parms[57] * x44 + x47 * x49;
	double x51 = parms[40] * x26 + parms[41] * x34 + parms[42] * x12 + parms[47] * x21 + parms[48] * x44 + x29 * x50 + x45 * x46;
	double x52 = parms[56] * x35 - parms[58] * x31 + parms[59] * x49;
	double x53 = -parms[57] * x35 + parms[58] * x30 + parms[59] * x40;
	double x54 = -x12;
	double x55 = parms[47] * x54 + parms[48] * x34 + parms[49] * x39 + x29 * x53 + x46 * x52;
	double x56 = -parms[57];
	double x57 = parms[52] * x31 + parms[54] * x30 + parms[55] * x35 + parms[56] * x49 + x40 * x56;
	double x58 = parms[41] * x26 + parms[43] * x34 + parms[44] * x12 + parms[46] * x22 + parms[48] * x39 - x57;
	double x59 = -x33;
	double x60 = parms[46] * x12 + parms[48] * x27 + parms[49] * x43 + parms[56] * x30 - parms[57] * x31 - parms[59] * x44;
	double x61 = -0.08 * x25;
	double x62 = parms[30] * x11 + parms[31] * x12 + parms[37] * x24 + parms[38] * x22 + x25 * x51 - 0.08 * x33 * x55 + x58 * x59 + x60 * x61;
	double x63 = parms[42] * x26 + parms[44] * x34 + parms[45] * x12 + parms[46] * x43 + parms[47] * x48 + x28 * x50 + x29 * x45;
	double x64 = parms[31] * x11 + parms[33] * x12 - parms[36] * x24 + parms[38] * x38 + x63;
	double x65 = -x60;
	double x66 = parms[36] * x54 + parms[37] * x11 + parms[39] * x24 + x25 * x65 + x55 * x59;
	double x67 = x25 * x55;
	double x68 = parms[38] * x12 + parms[39] * x38 + x33 * x65 + x67;
	double x69 = x10 * x68;
	double x70 = parms[38] * x32 + parms[39] * x21 + parms[46] * x35 + parms[47] * x26 + parms[49] * x21 + x28 * x53 + x29 * x52;
	double x71 = x70 * x8;
	double x72 = parms[21] * x9 + parms[23] * x6 - parms[26] * x23 + parms[28] * x19 + x10 * x64 + x37 * x62 - 0.17 * x66 + 0.005 * x69 + 0.005 * x71;
	double x73 = -x5;
	double x74 = x68 * x8;
	double x75 = x10 * x70;
	double x76 = parms[20] * x9 + parms[21] * x6 + parms[27] * x23 - parms[28] * x16 + x10 * x62 + x64 * x8 + 0.005 * x74 - 0.005 * x75;
	double x77 = parms[28] * x6 + parms[29] * x19 + x69 + x71;
	double x78 = x5 * x77;
	double x79 = parms[28] * x17 + parms[29] * x16 + x37 * x68 + x75;
	double x80 = x3 * x79;
	double x81 = 0.08 * x33;
	double x82 = -x25;
	double x83 = parms[32] * x11 + parms[34] * x12 + parms[36] * x21 + parms[37] * x42 + x51 * x59 + x58 * x82 + x60 * x81 - 0.08 * x67;
	double x84 = -parms[27];
	double x85 = parms[22] * x9 + parms[24] * x6 + parms[26] * x16 + x19 * x84 - 0.17 * x74 + 0.17 * x75 + x83;
	double x86 = parms[12] * x0 + parms[14] * x2 + parms[16] * x15 + parms[17] * x14 + 0.175 * x78 + 0.175 * x80 + x85;
	double x87 = -x82;
	double x88 = x29 * x59;
	double x89 = parms[56] * x87 - parms[58] * x88;
	double x90 = 0.175 * x3;
	double x91 = x90 + 0.17;
	double x92 = 0.175 * x5;
	double x93 = x10 * x91 + x8 * x92;
	double x94 = x10 * x92 + x37 * x91;
	double x95 = x94 - 0.08;
	double x96 = x25 * x95;
	double x97 = x29 * x93 + x46 * x96;
	double x98 = parms[59] * x97 + x89;
	double x99 = parms[46] * x87 + parms[47] * x59;
	double x100 = x46 * x59;
	double x101 = -parms[57] * x87 + parms[58] * x100;
	double x102 = x28 * x93 + x29 * x96;
	double x103 = parms[59] * x102 + x101;
	double x104 = parms[36] + parms[39] * x93 + parms[49] * x93 + x103 * x28 + x29 * x98 + x99;
	double x105 = x10 * x104;
	double x106 = x59 * x95;
	double x107 = -parms[48] * x59;
	double x108 = -x106;
	double x109 = -parms[56] * x100 + parms[57] * x88;
	double x110 = parms[49] * x106 - parms[59] * x108 + x107 - x109;
	double x111 = -parms[37];
	double x112 = parms[48] * x82;
	double x113 = x25 * (parms[49] * x96 + x103 * x29 + x112 + x46 * x98);
	double x114 = parms[39] * x94 + x110 * x59 + x111 + x113;
	double x115 = parms[41] * x59 + parms[43] * x82;
	double x116 = parms[52] * x88 + parms[54] * x100 + parms[55] * x87;
	double x117 = parms[56] * x97 + x102 * x56 + x116;
	double x118 = parms[50] * x88 + parms[51] * x100 + parms[52] * x87;
	double x119 = parms[57] * x108 + x118 + x47 * x97;
	double x120 = parms[40] * x59 + parms[41] * x82;
	double x121 = parms[51] * x88 + parms[53] * x100 + parms[54] * x87;
	double x122 = parms[58] * x102 + x108 * x41 + x121;
	double x123 = parms[35] + parms[36] * x93 + x110 * x81 + x111 * x94 - 0.08 * x113 + x59 * (parms[47] * x93 + parms[48] * x108 + x119 * x29 + x120 + x122 * x46) + x82 * (-parms[46] * x93 + parms[48] * x96 + x115 - x117);
	double x124 = parms[25] + parms[26] * x90 + 0.17 * x105 - 0.17 * x114 * x8 + x123 + x84 * x92;
	double x125 = -parms[47];
	double x126 = parms[42] * x59 + parms[44] * x82;
	double x127 = parms[46] * x106 + x119 * x28 + x122 * x29 + x125 * x96 + x126;
	double x128 = -0.17 * x8;
	double x129 = x128 - 0.08;
	double x130 = x129 * x59;
	double x131 = -x130;
	double x132 = parms[49] * x130 - parms[59] * x131 + x107 - x109;
	double x133 = 0.17 * x10;
	double x134 = x129 * x25;
	double x135 = x133 * x28 + x134 * x29;
	double x136 = parms[58] * x135 + x121 + x131 * x41;
	double x137 = x133 * x29 + x134 * x46;
	double x138 = parms[57] * x131 + x118 + x137 * x47;
	double x139 = parms[56] * x137 + x116 + x135 * x56;
	double x140 = parms[59] * x135 + x101;
	double x141 = parms[59] * x137 + x89;
	double x142 = parms[49] * x134 + x112 + x140 * x29 + x141 * x46;
	double x143 = parms[35] + parms[36] * x133 + x111 * x128 + x132 * x81 + x142 * x61 + x59 * (parms[47] * x133 + parms[48] * x131 + x120 + x136 * x46 + x138 * x29) + x82 * (-parms[46] * x133 + parms[48] * x134 + x115 - x139);
	double x144 = parms[46] * x130 + x125 * x134 + x126 + x136 * x29 + x138 * x28;
	double x145 = x29 * x61;
	double x146 = x46 * x61;
	double x147 = -x81;
	double x148 = parms[57] * x147 + x118 + x146 * x47;
	double x149 = parms[58] * x145 + x121 + x147 * x41;
	double x150 = parms[56] * x146 + x116 + x145 * x56;
	double x151 = parms[46] * x81 + x125 * x61 + x126 + x148 * x28 + x149 * x29;
	double x152 = parms[52] * x28 + parms[54] * x29;
	//
	M[0] = parms[5] + x0 * (parms[10] * x0 + parms[11] * x2 - parms[18] * x15 + x3 * x76 + x72 * x73) + x13 * (parms[18] * x2 + parms[19] * x13 + x3 * x77 + x73 * x79) + x15 * (parms[18] * x4 + parms[19] * x15 + x78 + x80) + x2 * (parms[11] * x0 + parms[13] * x2 + parms[18] * x13 - 0.175 * parms[26] * x7 - 0.175 * parms[27] * x9 - 0.175 * parms[29] * x23 + x3 * x72 + x5 * x76 - 0.175 * x66);
	M[1] = x86;
	M[2] = x85;
	M[3] = x83;
	M[4] = x63;
	M[5] = x57;
	M[6] = x86;
	M[7] = parms[15] + x124 + x90 * (parms[26] + parms[29] * x90 + x105 + x114 * x37) + x92 * (parms[29] * x92 + x10 * x114 + x104 * x8 + x84);
	M[8] = x124;
	M[9] = x123;
	M[10] = x127;
	M[11] = x117;
	M[12] = x85;
	M[13] = x124;
	M[14] = parms[25] + x128 * (parms[39] * x128 + x111 + x132 * x59 + x142 * x25) + x133 * (parms[36] + parms[39] * x133 + parms[49] * x133 + x140 * x28 + x141 * x29 + x99) + x143;
	M[15] = x143;
	M[16] = x144;
	M[17] = x139;
	M[18] = x83;
	M[19] = x123;
	M[20] = x143;
	M[21] = parms[35] + x59 * (parms[48] * x147 + x120 + x148 * x29 + x149 * x46) + x61 * (parms[49] * x61 + x112 + x29 * (parms[59] * x145 + x101) + x46 * (parms[59] * x146 + x89)) + x81 * (parms[49] * x81 - parms[59] * x147 + x107 - x109) + x82 * (parms[48] * x61 + x115 - x150);
	M[22] = x151;
	M[23] = x150;
	M[24] = x63;
	M[25] = x127;
	M[26] = x144;
	M[27] = x151;
	M[28] = parms[45] + x28 * (parms[50] * x28 + parms[51] * x29) + x29 * (parms[51] * x28 + parms[53] * x29);
	M[29] = x152;
	M[30] = x57;
	M[31] = x117;
	M[32] = x139;
	M[33] = x150;
	M[34] = x152;
	M[35] = parms[55];
	//
	return M;
}


Eigen::Matrix<double, DOF, DOF> sym_M_gluon(const double* parms, const double* q)
{
	Eigen::Matrix<double, DOF* DOF, 1> M_temp;
	Eigen::Matrix<double, DOF, DOF> M;
	M_temp = M_gluon(parms, q);
	for (int i = 0; i < DOF; ++i)
	{
		M.row(i) = M_temp.block(i * DOF, 0, DOF, 1).transpose();
		std::cout << M_temp.block(i * DOF, 0, DOF, 1) << std::endl;
	}
	return M;
}