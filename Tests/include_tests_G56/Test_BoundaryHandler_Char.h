// In the following routines we test the boundary handler for characteristic boundary conditions
namespace Test_BoundaryHanlder_Char
{
	using namespace dealii;
	#include "basic_routines.h"	

	//The matrix X_minus will not necessarily be equal to the one obtained from Mathematica since the eigenvectors
	//are unique only upto a multiplication with a scalar. So we check the matrix Xminus.(B.Xminus).inverse
	TEST(XminusBXminusInv,HandlesXminusBXminusInv)
	{
		const unsigned int dim = 2;

		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

		if (constants.constants.bc_type == characteristic)
		{
			MatrixOpt::Base_MatrixOpt matrix_opt;
			Full_matrix Xminus = matrix_opt.compute_Xminus(G26.system_data.Ax.matrix,
													  		G26.constants.nBC);

			Full_matrix XminusBXminusInv = Xminus * G26.B_tilde_inv;

			// result from Mathematica
			Full_matrix XminusBXminusInv_result;
			XminusBXminusInv_result.resize(G26.constants.nEqn,G26.constants.nBC);

			XminusBXminusInv_result << -0.7320760845083695,0.,-0.04258558465265678,0.17962335320178804,0.,0.,
   0.9999848117987198,0.,2.035415552714692e-6,-4.009815425254359e-6,0.,0.,0.,
   -0.5329013214962034,0.,0.,0.,-0.336674900590287,0.44635919602448737,0.,
   -0.2765787871923419,-0.363286456761382,0.,0.,-0.2986064462398187,0.,
   0.014355332434656191,-0.6202931925639745,0.,0.,0.,0.5166114566779891,0.,0.,0.,
   0.05105511219630321,0.1493032231199093,0.,-0.007177666217328071,-0.1423165533358946,
   -0.9049262992357638,0.,0.08585693972219449,0.,0.4790398184551206,
   0.054174968443899046,0.,0.,0.,0.03142090312343626,0.,0.,0.,-0.8089240169116636,
   -0.08342817269378511,0.,0.007828238411497532,0.6002454616153075,0.,0.,0.,
   -0.2579073094720947,0.,0.,0.,0.08161866336951966,0.041714086346892595,0.,
   -0.003914119205748783,-0.0039160690502147155,0.5924133235148781,0.,0.,
   0.19343048210407104,0.,0.,0.,-0.06121399752713974,0.13539697851655969,0.,
   0.27847970150368445,0.19375493517827924,0.,0.,-0.3027526010282279,0.,
   -0.22130361634744403,0.28837871164721884,0.,0.,0.,0.017776013478002783,0.,0.,0.,
   0.6274853895309347,0.15137630051411397,0.,0.110651808173722,0.09766237074935497,
   0.48370345314592883,0.;

			Compare_Float_Mat(XminusBXminusInv_result,XminusBXminusInv);			
		}

	}

	TEST(BhatCheck,HandlesBhat)
	{
		const unsigned int dim = 2;
		
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

		if (constants.constants.bc_type == characteristic)
		{
			Full_matrix B_hat;
			const unsigned int nEqn = G26.constants.nEqn;
			B_hat.resize(nEqn,nEqn);

			B_hat << 7.320760845083695e-6,-0.7320760845083695,0.,0.014159084261081336,
   					-0.11611968880077986,0.,0.,-0.03774049175533217,0.,0.15918705204753542,0.,0.,0.,
   					-0.012669616290574039,0.07518841723210296,0.,0.,-9.999848117987199e-6,
   					0.9999848117987198,0.,7.056433472442031e-6,-0.000011206465507567663,0.,0.,
   1.8038400673010677e-6,0.,-3.553606395956214e-6,0.,0.,0.,9.913905452287864e-7,
   -2.5952535768521623e-6,0.,0.,0.,0.,0.31144073677435274,0.,0.,-0.47227149971929205,
   0.,0.,0.0722391263797172,0.,0.46497036917737555,0.,0.,0.,0.,-0.29837036202724293,0.,
   -4.463591960244874e-6,0.44635919602448737,0.,0.23623920768302387,0.1201673004163438,
   0.,0.,-0.2451115682189218,0.,-0.32195423963426584,0.,0.,0.,-0.21129551118530066,
   0.15447613461748289,0.,0.,2.986064462398187e-6,-0.2986064462398187,0.,
   0.04799571240009226,0.35907086726996135,0.,0.,0.01272208212741761,0.,
   -0.5497205289252467,0.,0.,0.,-0.04293085093520665,-0.14749154209425064,0.,0.,0.,0.,
   -0.26512825468477563,0.,0.,0.457834982905405,0.,0.,0.057951248507901384,0.,
   -0.5280047346554112,0.,0.,0.,0.,0.045246415110375394,0.,-1.493032231199093e-6,
   0.1493032231199093,0.,-0.023997856200046144,0.00518187369441464,0.,
   0.3694346146587908,-0.0063610410637087835,0.,-0.1261247615039033,0.,
   -0.8019700519330534,0.,0.021465425467603345,-0.02498978287202104,0.,
   -0.19747110783829272,-8.58569397221945e-7,0.08585693972219449,0.,
   -0.3549481833367213,0.1205209038682709,0.,0.,0.42453798547892885,0.,
   0.04801131572053684,0.,0.,0.,0.3174759335828215,-0.3930405934342922,0.,0.,0.,0.,
   0.09238656712458097,0.,0.,0.027846050370031877,0.,0.,0.3809841470552141,0.,
   -0.25994153867353,0.,0.,0.,0.,-0.7168902444325184,0.,8.342817269378512e-7,
   -0.08342817269378511,0.,-0.06230930304413617,-0.34059217079115,0.,0.,
   0.0069375956591326395,0.,0.5319536899643054,0.,0.,0.,0.05573052562170101,
   0.12436822925691537,0.,0.,0.,0.,0.1180469020683624,0.,0.,-0.22856440192524746,0.,0.,
   -0.07871792449348913,0.,0.2936456857345018,0.,0.,0.,0.,0.07233265709753305,0.,
   -4.17140863468926e-7,0.041714086346892595,0.,0.031154651522068096,
   0.04937022210909286,0.,-0.24185172657296442,-0.003468797829566335,0.,
   -0.0034705258342324894,0.,0.5250126382958404,0.,-0.027865262810850524,
   0.0024534781740968644,0.,0.12927518560510912,0.,0.,-0.0885351765512718,0.,0.,
   0.1714233014439356,0.,0.,0.059038443370116836,0.,-0.22023426430087634,0.,0.,0.,0.,
   -0.05424949282314978,0.,-1.353969785165597e-6,0.13539697851655969,0.,
   -0.2216391040565203,-0.022679070392543937,0.,0.,0.24679620966461205,0.,
   0.17171084049434482,0.,0.,0.,0.1982410300567362,-0.19307751140395063,0.,0.,
   3.027526010282279e-6,-0.3027526010282279,0.,0.13442623266923762,
   -0.23480022909126316,0.,0.,-0.19612522350717201,0.,0.2555689789891422,0.,0.,0.,
   -0.12023668867543305,0.2499650869874793,0.,0.,0.,0.,-0.09273926898297088,0.,0.,
   0.015753581771417193,0.,0.,-0.28886643619271046,0.,0.1573900503814748,0.,0.,0.,0.,
   0.5560944475305264,0.,-1.5137630051411398e-6,0.15137630051411397,0.,
   -0.06721311633461881,0.018664560626485263,0.,-0.19747110783829275,0.098062611753586,
   0.,0.08655102256162821,0.,0.42867102411239866,0.,0.06011834433771651,
   -0.07220616997105525,0.,0.10555274704536878;


			Compare_Float_Mat(B_hat,G26.B_hat);			
		}
	}

	TEST(BCrhs,HandlesBCrhsSystemA)
	{
		const unsigned int dim = 2;
		std::string folder_name = "../system_matrices/";
		Constants::Base_Constants constants(input_file);
		G26::G26<dim> G26(constants.constants,folder_name);

 
		const double  theta0 = 1.0;
		const double  theta1 = 2.0;
		const double  theta2 = 3.0;
		const double  theta3 = 4.0;

 
		const double  vn0 = 1.0;
		const double  vn1 = 2.0;
		const double  vn2 = 3.0;
		const double  vn3 = 4.0;

 
		const double  vt0 = 1.0;
		const double  vt1 = 2.0;
		const double  vt2 = 3.0;
		const double  vt3 = 4.0;

		
		if (constants.constants.bc_type == characteristic)
		{
			// for this particular system, p and normal_vector do not matter. Everything regulated by b_id
			// consider the wall with b_id == 0
			Tensor<1,dim> p;
			Tensor<1,dim> normal_vector;
			unsigned int b_id = 0;
			Vector<double> bc_rhs_manuel(G26.constants.nBC);
			Vector<double> bc_rhs(G26.constants.nBC);
			bc_rhs_manuel = 0;

			
			bc_rhs_manuel(0) += vn0;
			bc_rhs_manuel(1) += -0.5 * vt0;
			bc_rhs_manuel(2) += -2*sqrt(2/15.0) * (-sqrt(3/2.0)) * theta0;
			bc_rhs_manuel(3) += -sqrt(2)/15.0 * (-sqrt(3/2.0)) * theta0;
			bc_rhs_manuel(4) += 1/(15.0 * sqrt(2)) * (-sqrt(3/2.0)) * theta0;
			bc_rhs_manuel(5) += -1/(2*sqrt(14)) * vt0;

			G26.build_BCrhs(p,normal_vector,bc_rhs,b_id);

			
			for (int i = 0 ; i < G26.constants.nBC ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "Failed at" << i;


			
			b_id = 1;
			bc_rhs_manuel = 0;
			bc_rhs = 0;

			bc_rhs_manuel(0) += vn1;
			bc_rhs_manuel(1) += -0.5 * vt1;
			bc_rhs_manuel(2) += -2*sqrt(2/15.0) * (-sqrt(3/2.0)) * theta1;
			bc_rhs_manuel(3) += -sqrt(2)/15.0 * (-sqrt(3/2.0)) * theta1;
			bc_rhs_manuel(4) += 1/(15.0 * sqrt(2)) * (-sqrt(3/2.0)) * theta1;
			bc_rhs_manuel(5) += -1/(2*sqrt(14)) * vt1;

			G26.build_BCrhs(p,normal_vector,bc_rhs,b_id);

			
			for (int i = 0 ; i < G26.constants.nBC ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "Failed at" << i;

			b_id = 2;
			bc_rhs_manuel = 0;
			bc_rhs = 0;

			
			bc_rhs_manuel(0) += vn2;
			bc_rhs_manuel(1) += -0.5 * vt2;
			bc_rhs_manuel(2) += -2*sqrt(2/15.0) * (-sqrt(3/2.0)) * theta2;
			bc_rhs_manuel(3) += -sqrt(2)/15.0 * (-sqrt(3/2.0)) * theta2;
			bc_rhs_manuel(4) += 1/(15.0 * sqrt(2)) * (-sqrt(3/2.0)) * theta2;
			bc_rhs_manuel(5) += -1/(2*sqrt(14)) * vt2;

			G26.build_BCrhs(p,normal_vector,bc_rhs,b_id);

			
			for (int i = 0 ; i < G26.constants.nBC ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "Failed at" << i;


			b_id = 3;
			bc_rhs_manuel = 0;
			bc_rhs = 0;

			

			bc_rhs_manuel(0) += vn3;
			bc_rhs_manuel(1) += -0.5 * vt3;
			bc_rhs_manuel(2) += -2*sqrt(2/15.0) * (-sqrt(3/2.0)) * theta3;
			bc_rhs_manuel(3) += -sqrt(2)/15.0 * (-sqrt(3/2.0)) * theta3;
			bc_rhs_manuel(4) += 1/(15.0 * sqrt(2)) * (-sqrt(3/2.0)) * theta3;
			bc_rhs_manuel(5) += -1/(2*sqrt(14)) * vt3;

			G26.build_BCrhs(p,normal_vector,bc_rhs,b_id);

			for (int i = 0 ; i < G26.constants.nBC ; i++)
				EXPECT_NEAR(bc_rhs(i),bc_rhs_manuel(i),1e-5) << "Failed at" << i;
		}
		
	}
}