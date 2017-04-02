// exact solution for poisson heat conduction
namespace ExactSolution
{
	using namespace dealii;

	template<int dim>
	class
	PoissonHeat:public Base_ExactSolution<dim>
	{
	public:
		PoissonHeat(const constant_numerics &constants,const Sparse_matrix &S_half,const int nEqn,const int Ntensors);

		virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

		virtual double s_r(const double r,const double phi) const;
		virtual double s_phi(const double r,const double phi) const;
		virtual double thetaP(const double r,const double phi) const;	
		virtual double R_rr(const double r,const double phi) const;
		virtual double R_thetatheta(const double r,const double phi) const;
		virtual double R_rtheta(const double r,const double phi) const;
		virtual double R_zz(const double r,const double phi) const;
	};

	template<int dim>
	PoissonHeat<dim>::PoissonHeat(const constant_numerics &constants,const Sparse_matrix &S_half,const int nEqn,
								const int Ntensors)
	:
	Base_ExactSolution<dim>(constants,S_half,nEqn,Ntensors)
	{

	}


	template<>
	void
	PoissonHeat<2>::vector_value(const Point<2> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha - 0.816496580927726) < 1e-5,ExcMessage("Exact Solution does not correspond to the given value of alpha"));

		// variables for which we need the 
		const unsigned int ID_theta = this->constants.variable_map.find("theta")->second;
		const unsigned int ID_heat = this->constants.variable_map.find("qy")->second;
		const unsigned int ID_stress = this->constants.variable_map.find("sigmayy")->second;
		const double y = p(1);
		bool developed_exact_solution = false;

		value = 0;

		Assert(ID_theta == 3,ExcMessage("Wrong ID for theta"));

		// G20
		if (this->Ntensors == 6)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.2829983743239428 + 0.0003362282005094463*cosh(7.453559924999298*y) - 
					0.026127890589687234*pow(y,2) + 0.408248290463863*pow(y,4);

					value[ID_stress] = -(pow(5,0.5)*(-0.0013576450198781716 + 
						0.0004853036051831808*cosh(7.453559924999298*y) - 0.03771236166328255*pow(y,2)
						))/2.;


				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.2898753498434943 + 0.02291602260455048*cosh(2.484519974999766*y) - 
					0.0783836717690617*pow(y,2) + 0.13608276348795437*pow(y,4);

					value[ID_stress] = -(pow(5,0.5)*(-0.03665641553671062 + 0.03307642954873192*cosh(2.484519974999766*y) - 
						0.1131370849898476*pow(y,2)))/2.;

				}			
			}
		}

		// G26
		if (this->Ntensors == 8)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{

					developed_exact_solution = true;

					value[ID_theta] = -1.2852476333586613 - 0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) + 
					0.0015643311536178803*cosh(6.573421981221795 * y) - 
					2.168404344971009e-19*sinh(6.573421981221795 * y);

					value[ID_stress] = -0.002715290039756343 - 0.03771236166328255*pow(y,2) + 
					0.0011289587658037511*cosh(6.573421981221795*y);

					value[ID_heat] = -0.21081851067789195*pow(y,3);

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.3650564316627332 - 0.548685702383432*pow(y,2) + 0.13608276348795437*pow(y,4) + 
					0.09338201417693054*cosh(2.191140660407265*y);

					value[ID_stress] = -0.07331283107342124 - 0.1131370849898476*pow(y,2) + 
					0.06739266377815037*cosh(2.191140660407265*y);

					value[ID_heat] = -0.21081851067789195*pow(y,3);
				}

				if (fabs(this->constants.tau - 0.5) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7298151776227189 - 0.9144761706390533*pow(y,2) + 0.08164965809277261*pow(y,4) + 
					0.45942750968124524*cosh(1.314684396244359*y) - 
					2.7755575615628914e-17*sinh(1.314684396244359*y);

					value[ID_stress] = -0.3394112549695428 - 0.1885618083164127*pow(y,2) + 
					0.33156324548448296*cosh(1.314684396244359*y);

					value[ID_heat] = -0.21081851067789195*pow(y,3);
				}			
			}
		}

		// G35
		if (this->Ntensors == 9)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.2857735434320117 + 0.002601070885430934*cosh(4.662524041201569*y) - 
					0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) + 
					2.168404344971009e-19*sinh(4.662524041201569*y);

					value[ID_stress] = -(pow(5,0.5)*(-0.002715290039756343 + 0.00187716121985606*cosh(4.662524041201569*y) - 
						0.03771236166328255*pow(y,2) - 
						1.0842021724855044e-19*sinh(4.662524041201569*y)))/2.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.3660400339849652 + 0.09916339949943631*cosh(1.554174680400523*y) - 
					0.548685702383432*pow(y,2) + 0.13608276348795437*pow(y,4) + 
					6.938893903907228e-18*sinh(1.554174680400523*y);

					value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.07156501924344744*cosh(1.554174680400523*y) - 
						0.1131370849898476*pow(y,2)))/2.;

				}		
			}
		}

		// G45
		if (this->Ntensors == 11)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.2995436998701937 + 0.012092296117067552*cosh(3.849386397950811*y) + 
					0.0015393861867350283*cosh(5.674996469312898*y) + 
					0.000022293298256729303*cosh(10.161219234535377*y) - 
					0.18289523412781067*pow(y,2) + 0.408248290463863*pow(y,4) - 
					8.673617379884035e-19*sinh(3.849386397950811*y) - 
					5.421010862427522e-19*sinh(5.674996469312898*y) - 
					8.470329472543003e-21*sinh(10.161219234535377*y);

					value[ID_stress] = -(pow(5,0.5)*(-0.002715290039756343 + 
						0.0031225329963700333*cosh(3.849386397950811*y) - 
						0.0007338661922650514*cosh(5.674996469312898*y) + 
						0.000033602606211988805*cosh(10.161219234535377*y) - 
						0.03771236166328255*pow(y,2) - 
						1.5178830414797062e-18*sinh(3.849386397950811*y) - 
						1.0842021724855044e-19*sinh(5.674996469312898*y) + 
						1.0164395367051604e-20*sinh(10.161219234535377*y)))/2.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7443580752405021 + 0.39756135791713554*cosh(1.283128799316937*y) + 
					0.07237762337695752*cosh(1.8916654897709662*y) + 
					0.003071248679976172*cosh(3.3870730781784593*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4) - 2.7755575615628914e-17*sinh(1.283128799316937*y) - 
					3.469446951953614e-17*sinh(1.8916654897709662*y) - 
					1.5178830414797062e-18*sinh(3.3870730781784593*y);

					value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.10266027610966066*cosh(1.283128799316937*y) - 
						0.03450433122665434*cosh(1.8916654897709662*y) + 
						0.004629281804058666*cosh(3.3870730781784593*y) - 
						0.1131370849898476*pow(y,2) - 
						4.85722573273506e-17*sinh(1.283128799316937*y) - 
						6.938893903907228e-18*sinh(1.8916654897709662*y) + 
						8.673617379884035e-19*sinh(3.3870730781784593*y)))/2.;

				}		
			}
		}

		// G56
		if (this->Ntensors == 12)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.2991542242917795 + 0.008670378733644997*cosh(3.536982581771994*y) + 
					0.004430725803158755*cosh(5.204696760361044*y) + 
					0.000056930062888828204*cosh(10.45417499077846*y) - 0.18289523412781067*pow(y,2) + 
					0.408248290463863*pow(y,4) - 1.3010426069826053e-18*sinh(5.204696760361044*y) + 
					1.0164395367051604e-20*sinh(10.45417499077846*y);

					value[ID_stress] = -0.002715290039756343 + 0.0026403019195414515*cosh(3.536982581771994*y) - 
					0.0006407561824316862*cosh(5.204696760361044*y) + 
					0.00008292874116147389*cosh(10.45417499077846*y) - 0.03771236166328255*pow(y,2) - 
					3.469446951953614e-18*sinh(3.536982581771994*y) - 
					3.2526065174565133e-18*sinh(5.204696760361044*y);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7435149236809255 + 0.27723853580237484*cosh(1.178994193923998*y) + 
					0.18476156620816786*cosh(1.734898920120348*y) + 
					0.008511751878547334*cosh(3.4847249969261536*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4) - 5.551115123125783e-17*sinh(1.734898920120348*y) + 
					8.673617379884035e-19*sinh(3.4847249969261536*y);

					value[ID_stress] = -(pow(5,0.5)*(-0.07331283107342124 + 0.0844246209694861*cosh(1.178994193923998*y) - 
						0.02671957622366166*cosh(1.734898920120348*y) + 
						0.0123988773689772*cosh(3.4847249969261536*y) - 0.1131370849898476*pow(y,2) - 
						1.1102230246251565e-16*sinh(1.178994193923998*y) - 
						1.3357370765021415e-16*sinh(1.734898920120348*y)))/2.;

				}		
			}
		}

		// G71
		if (this->Ntensors == 15)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.3012568770728508 + 0.014368644008044136*cosh(2.9184816217908978*y) + 
					0.000847741867453309*cosh(3.885676724821849*y) + 
					0.0004461674640337943*cosh(5.896124120768356*y) + 
					0.00014238973643247108*cosh(9.591077444821629*y) - 0.18289523412781067*pow(y,2) + 
					0.408248290463863*pow(y,4) + 8.673617379884035e-19*sinh(3.885676724821849*y) - 
					1.3552527156068805e-20*sinh(9.591077444821629*y);

					value[ID_stress] = -0.002715290039756343 + 0.0018696658823278787*cosh(2.9184816217908978*y) - 
					0.0006185894211626305*cosh(3.885676724821849*y) + 
					0.0004424972314741653*cosh(5.896124120768356*y) + 
					0.00010260159078614521*cosh(9.591077444821629*y) - 0.03771236166328255*pow(y,2) - 
					4.445228907190568e-18*sinh(2.9184816217908978*y) + 
					3.7947076036992655e-19*sinh(3.885676724821849*y) + 
					2.439454888092385e-19*sinh(5.896124120768356*y) - 
					4.0657581468206416e-20*sinh(9.591077444821629*y);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.747245655093501 + 0.41369103486028996*cosh(0.9728272072636326*y) + 
					0.02470452567666627*cosh(1.2952255749406163*y) + 
					0.021023428196056175*cosh(1.9653747069227854*y) + 
					0.017552620592354434*cosh(3.197025814940543*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4) + 2.6020852139652106e-17*sinh(1.2952255749406163*y) - 
					1.734723475976807e-18*sinh(3.197025814940543*y);

					value[ID_stress] = -0.07331283107342124 + 0.053829993510186594*cosh(0.9728272072636326*y) - 
					0.01802666451326118*cosh(1.2952255749406163*y) + 
					0.020850486695610182*cosh(1.9653747069227854*y) + 
					0.012647869434713875*cosh(3.197025814940543*y) - 0.1131370849898476*pow(y,2) - 
					1.249000902703301e-16*sinh(0.9728272072636326*y) + 
					1.214306433183765e-17*sinh(1.2952255749406163*y) + 
					1.214306433183765e-17*sinh(1.9653747069227854*y) - 
					6.071532165918825e-18*sinh(3.197025814940543*y);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}		
			}
		}

		// G84
		if (this->Ntensors == 16)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.3016616131237297 + 0.012498585238508764*cosh(2.867124849692578*y) + 
					0.0026996208249191254*cosh(3.6131625584441966*y) + 
					0.0008429866277644128*cosh(6.002769036811866*y) - 0.18289523412781067*pow(y,2) + 
					0.408248290463863*pow(y,4) - 2.6020852139652106e-18*sinh(2.867124849692578*y) - 
					8.673617379884035e-19*sinh(3.6131625584441966*y) + 
					1.0842021724855044e-19*sinh(6.002769036811866*y);

					value[ID_stress] = -0.002715290039756343 + 0.0018033048153696892*cosh(2.867124849692578*y) - 
					0.0007033685726265145*cosh(3.6131625584441966*y) + 
					0.0008142812934811178*cosh(6.002769036811866*y) - 0.03771236166328255*pow(y,2) + 
					1.951563910473908e-18*sinh(2.867124849692578*y) - 
					1.6263032587282567e-18*sinh(3.6131625584441966*y) + 
					1.0842021724855044e-19*sinh(6.002769036811866*y);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7479896311184409 + 0.36064751628260816*cosh(0.9557082832308594*y) + 
					0.07833003503050456*cosh(1.2043875194813989*y) + 
					0.04090384023127845*cosh(2.000923012270622*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4) - 8.326672684688674e-17*sinh(0.9557082832308594*y) - 
					3.469446951953614e-17*sinh(1.2043875194813989*y) + 
					3.469446951953614e-18*sinh(2.000923012270622*y);

					value[ID_stress] = -0.07331283107342124 + 0.05203448153153864*cosh(0.9557082832308594*y) - 
					0.0204083789933134*cosh(1.2043875194813989*y) + 
					0.0395109849134863*cosh(2.000923012270622*y) - 0.1131370849898476*pow(y,2) + 
					5.551115123125783e-17*sinh(0.9557082832308594*y) - 
					4.683753385137379e-17*sinh(1.2043875194813989*y) + 
					6.938893903907228e-18*sinh(2.000923012270622*y);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}		
			}
		}

		// G105
		if (this->Ntensors == 19)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{

					developed_exact_solution = true;

					value[ID_theta] = -1.2997563950037463 + 0.009276947895287067*cosh(2.492986999351348*y) + 
					0.00034793324875276896*cosh(3.170853933340145*y) + 
					0.0019572547988298645*cosh(4.388700162109401*y) + 
					0.0017525338749780097*cosh(5.365272100745118*y) + 
					0.0006825478821284116*cosh(7.820516007402679*y) + 
					6.90472273800691e-7*cosh(10.44464164243648*y) + 
					7.880009718028203e-6*cosh(13.163405846932237*y) - 0.18289523412781067*pow(y,2) + 
					0.408248290463863*pow(y,4);

					value[ID_stress] = -0.002715290039756343 + 0.0008236716794642021*cosh(2.492986999351348*y) - 
					0.00023601868281006479*cosh(3.170853933340145*y) + 
					0.0011001213990571828*cosh(4.388700162109401*y) + 
					0.0006230403314441149*cosh(5.365272100745118*y) - 
					0.00020406994398414134*cosh(7.820516007402679*y) + 
					2.59443373639468e-6*cosh(10.44464164243648*y) + 
					0.000010993753807602042*cosh(13.163405846932237*y) - 0.03771236166328255*pow(y,2);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7441170080467425 + 0.2660463749885166*cosh(0.8309956664504494*y) + 
					0.00981560623079519*cosh(1.0569513111133817*y) + 
					0.07058859198806446*cosh(1.462900054036467*y) + 
					0.07231192665512148*cosh(1.7884240335817059*y) + 
					0.049617730356627546*cosh(2.6068386691342265*y) + 
					0.00009537413540004746*cosh(3.4815472141454937*y) + 
					0.0025738046900230706*cosh(4.387801948977413*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4);

					value[ID_stress] = -0.07331283107342124 + 0.023621439613073676*cosh(0.8309956664504494*y) - 
					0.006658364677359996*cosh(1.0569513111133817*y) + 
					0.03967598936112579*cosh(1.462900054036467*y) + 
					0.025707489820209463*cosh(1.7884240335817059*y) - 
					0.014834838286982257*cosh(2.6068386691342265*y) + 
					0.0003583661268530281*cosh(3.4815472141454937*y) + 
					0.0035908299765454493*cosh(4.387801948977413*y) - 0.1131370849898476*pow(y,2);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}		
			}
		}

		// G120
		if (this->Ntensors == 20)
		{
			if(fabs(this->constants.alpha-0.816496580927726) < 1e-5)
			{
				if (fabs(this->constants.tau - 0.1) < 1e-5 )
				{

					developed_exact_solution = true;

					value[ID_theta] = -1.2995390567873166 + 0.008828936885305286*cosh(2.480359051432615*y) + 
					0.0008633188573549727*cosh(3.028788951922567*y) + 
					0.003055077322774268*cosh(4.434798678173429*y) + 
					0.0011530952298520565*cosh(7.524303732560254*y) + 
					0.000013607923731562518*cosh(12.976581617858354*y) - 0.18289523412781067*pow(y,2) + 
					0.408248290463863*pow(y,4);

					value[ID_stress] = -0.002715290039756343 + 0.0008187280616851396*cosh(2.480359051432615*y) - 
					0.0002604386379633425*cosh(3.028788951922567*y) + 
					0.0016721283891165937*cosh(4.434798678173429*y) - 
					0.00017655344915385124*cosh(7.524303732560254*y) + 
					0.00001945869320751019*cosh(12.976581617858354*y) - 0.03771236166328255*pow(y,2);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}

				if (fabs(this->constants.tau - 0.3) < 1e-5)
				{
					developed_exact_solution = true;

					value[ID_theta] = -1.7436949415826675 + 0.2523069701617391*cosh(0.8267863504775383*y) + 
					0.02533033765112578*cosh(1.0095963173075224*y) + 
					0.11056926376618245*cosh(1.4782662260578097*y) + 
					0.07739730163732568*cosh(2.508101244186751*y) + 
					0.004179811421577409*cosh(4.325527205952785*y) - 0.548685702383432*pow(y,2) + 
					0.13608276348795437*pow(y,4);

					value[ID_stress] = -0.07331283107342124 + 0.02339701815900208*cosh(0.8267863504775383*y) - 
					0.007641439290718827*cosh(1.0095963173075224*y) + 
					0.060517618827160266*cosh(1.4782662260578097*y) - 
					0.011850504802646761*cosh(2.508101244186751*y) + 
					0.005976934448057989*cosh(4.325527205952785*y) - 0.1131370849898476*pow(y,2);

					value[ID_heat] = -(pow(0.4,0.5)*pow(y,3))/3.;

				}		
			}
		}


		Assert(developed_exact_solution,ExcMessage("Exact solution not initialized"));

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
		MatrixOpt::Base_MatrixOpt matrix_opt;
		value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	template<>
	void
	PoissonHeat<1>::vector_value(const Point<1> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha) < 1e-5,ExcMessage("Exact Solution does not correspond to the given value of alpha"));
		Assert(fabs(this->constants.theta0+1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.theta1-1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.tau-0.1) < 1e-5,ExcMessage("Incorrect tau value"));
		bool developed_exact_solution = false;

		// variables for which we need the exact solution
		const unsigned int ID_theta = this->constants.variable_map_1D.find("theta")->second;

		AssertDimension(ID_theta,2);
		const double x = p(0);

		value = 0;

		if(fabs(this->constants.alpha) < 1e-5)
			if (fabs(this->constants.tau - 0.1) < 1e-5 )
			{
				if (this->Ntensors == 6)
				{
					developed_exact_solution = true;
					value[ID_theta] = 0.9950211932019228*x + 0.00974532392134874*sinh(4.47213595499958*x);
				}
						

				if (this->Ntensors == 8)
				{
					developed_exact_solution = true;

					value[ID_theta] = 0.9918095457622201*x + 0.004077339347408366*sinh(2.517180972479634*x) + 
   								  0.0031858247584759585*sinh(6.71508535912671*x);	
				}
			

				if (this->Ntensors == 9)
				{
					developed_exact_solution = true;
					
					value[ID_theta] = 0.9861871317460368*x + 0.0009740846199919195*sinh(2.2483363346888074*x) + 
   									0.008298820931674668*sinh(4.010385905150765*x);					
				}


				if (this->Ntensors == 11)
				{

					developed_exact_solution = true;

					value[ID_theta] = 0.9857643470749898*x + 0.0004959003321893117*sinh(1.9303025585770384*x) + 
   								  0.002543347245910525*sinh(2.6658847234287704*x) 
   								  + 0.005290233872198295*sinh(4.9438923511630914*x);

				}

				if (this->Ntensors == 12)
				{

					developed_exact_solution = true;

					value[ID_theta] = 0.9899466270992661*x + 0.00010116292970230688*sinh(1.8186266940442457*x) + 
   									  0.0019797047973874405*sinh(2.3891354429537364*x) + 0.004978603330380985*sinh(3.687003638620778*x) + 
   									  0.00038639730251888393*sinh(10.604506967165621*x);

				}

				if (this->Ntensors == 14)
				{
					developed_exact_solution = true;
					value[ID_theta] = 0.989552444757239*x + 7.554239778848198e-6*sinh(1.64277980395364*x) + 
   			0.0004996637794229688*sinh(2.0483113377580193*x) + 0.0028382726390473777*sinh(2.6857066212139475*x) + 
   			0.004092620421355574*sinh(4.230009784703013*x) + 0.00015174242599143398*sinh(12.325411809002668*x);
				}

				if (this->Ntensors == 16)
				{
					developed_exact_solution = true;

					value[ID_theta] = 0.9892927396698521*x + 7.949524702960454e-7*sinh(1.5080981526429016*x) + 
   0.00007599110607032715*sinh(1.8280360908575142*x) + 0.0010765095377805595*sinh(2.2459528687383394*x) + 
   0.003195836632356749*sinh(2.9712749000654357*x) + 0.0032524332462132983*sinh(4.745522597484798*x) + 
   0.0000635644033471918*sinh(13.934429182642614*x);
				}

				if (this->Ntensors == 18)
				{
					developed_exact_solution = true;

					value[ID_theta] = 
0.9891103738688487*x - 3.4684362295131658e-9*sinh(1.4006653082978795*x) + 
   9.938811076145555e-6*sinh(1.664568126811443*x) + 0.00025689789010162757*sinh(1.9815191312437987*x) + 
   0.001615975388904252*sinh(2.431133592076258*x) + 0.0032265940137877307*sinh(3.2473676118677286*x) + 
   0.002561075578033159*sinh(5.235889289787754*x) + 0.000028099072925427022*sinh(15.450468880651211*x);
				}

				if (this->Ntensors == 20)
				{
					developed_exact_solution = true;

					value[ID_theta] = 
0.9889757964874321*x + 8.548836749354673e-9*sinh(1.312499845273545*x) + 
   6.425109739508264e-7*sinh(1.5359624765465445*x) + 0.000046769982279263285*sinh(1.7928504374446674*x) + 
   0.0005340677389950208*sinh(2.1202457820361236*x) + 0.002012515257218323*sinh(2.6095498761963127*x) + 
   0.003085495180863663*sinh(3.5141607592624338*x) + 0.0020173922087678656*sinh(5.703742588502402*x) + 
   0.000012998563891241376*sinh(16.88777884716506*x);
				}

				if (this->Ntensors == 22)
				{
					developed_exact_solution = true;

					value[ID_theta] = 

					0.9888725604237197*x - 1.1267507505527853e-9*sinh(1.238538873534564*x) + 
   1.0425594006751506e-7*sinh(1.4314385108767749*x) + 5.310480436339387e-6*sinh(1.6466429341157962*x) + 
   0.000130574387363048*sinh(1.9067974753769268*x) + 0.0008484475423778083*sinh(2.2508407141902804*x) + 
   0.002257045383117194*sinh(2.7829365376803588*x) + 0.0028667510847416398*sinh(3.771937211055834*x) + 
   0.0015959510715178827*sinh(6.151532393786605*x) + 6.251452487753951e-6*sinh(18.257360408445823*x);

				}

				if (this->Ntensors == 24)
				{
					developed_exact_solution = true;

					value[ID_theta] = 
0.9887909148897903*x + 2.1248298302441273e-10*sinh(1.1753915070425638*x) - 
   4.753332814364334e-9*sinh(1.3444672058183817*x) + 7.269127503526152e-7*sinh(1.5286628972393546*x) + 
   0.000021091722168993028*sinh(1.7445093738681328*x) + 0.00026554174814196423*sinh(2.0119506654634387*x) + 
   0.0011453546231374565*sinh(2.376383830762091*x) + 0.0023788884543610805*sinh(2.951828444318971*x) + 
   0.0026215690952677067*sinh(4.0211509594000825*x) + 0.0012699863181260631*sinh(6.581408804336385*x) + 
   3.1095921280567805e-6*sinh(19.56789817463244*x);
				}


				if (
					this->Ntensors == 15 ||
				    this->Ntensors == 19 )
				{
					AssertThrow(1 == 0,ExcMessage("Should not have reached here"));
				}

			}

			
			

			Assert(developed_exact_solution,ExcMessage("exact solution not created"));

		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
			MatrixOpt::Base_MatrixOpt matrix_opt;
			value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
		}
	// need to fix the following routines
	template<int dim>
		double
		PoissonHeat<dim>::s_r(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::s_phi(const double r,const double phi)const
		{

			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::thetaP(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::R_rr(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::R_thetatheta(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::R_rtheta(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat<dim>::R_zz(const double r,const double phi)const
		{
			return 0;

		}

	// exact solution obtained from a very high order moment method Ntensors = 30
	template<int dim>
		class
		PoissonHeat_HigherOrder:public PoissonHeat<dim>
		{
		public:

			PoissonHeat_HigherOrder(const constant_numerics &constants,const Sparse_matrix &S_half,const int nEqn,const int Ntensors);

			virtual void vector_value(const Point<dim> &p,Vector<double> &value) const ;

			virtual double s_r(const double r,const double phi) const;
			virtual double s_phi(const double r,const double phi) const;
			virtual double thetaP(const double r,const double phi) const;	
			virtual double R_rr(const double r,const double phi) const;
			virtual double R_thetatheta(const double r,const double phi) const;
			virtual double R_rtheta(const double r,const double phi) const;
			virtual double R_zz(const double r,const double phi) const;

		};

	template<int dim>
	PoissonHeat_HigherOrder<dim>::PoissonHeat_HigherOrder(const constant_numerics &constants,const Sparse_matrix &S_half,const int nEqn,
								const int Ntensors)
	:
	PoissonHeat<dim>(constants,S_half,nEqn,Ntensors)
	{

	}

	template<>
	void 
	PoissonHeat_HigherOrder<1>::vector_value(const Point<1> &p,Vector<double> &value) const
	{
		// first we check the size of the value vector
		Assert((int)value.size() == this->nEqn,ExcNotInitialized());
		Assert(fabs(this->constants.alpha) < 1e-5,ExcMessage("Exact Solution does not correspond to the given value of alpha"));
		Assert(fabs(this->constants.theta0+1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.theta1-1) < 1e-5,ExcMessage("Incorrect temperature value"));
		Assert(fabs(this->constants.tau-0.1) < 1e-5,ExcMessage("Incorrect tau value"));
		bool developed_exact_solution = false;

		// variables for which we need the exact solution
		const unsigned int ID_theta = this->constants.variable_map_1D.find("theta")->second;

		AssertDimension(ID_theta,2);
		const double x = p(0);

		value = 0;
		value[ID_theta] = 0.9886241450400803*x - 1.4284622829507996e-9*sinh(1.2779468215207301*x) + 
   						5.266949286451139e-8*sinh(1.4173548123589153*x) + 1.1902338825856906e-6*sinh(1.5772937256341433*x) + 
   						0.00002541376197682671*sinh(1.766211842724685*x) + 0.00020307410067078814*sinh(1.9973620316950804*x) + 
   0.0008206049409972356*sinh(2.29896310304997*x) + 0.0017197261477691537*sinh(2.7338379482754673*x) + 
   0.002324240918088127*sinh(3.4334119104667082*x) + 0.0019292768842024432*sinh(4.722835521599962*x) + 
   0.0006649745620110362*sinh(7.7807719198946925*x) + 4.5110234469975064e-7*sinh(23.209091470407046*x);

   		// The above values correspond to a unsymmetric system, therefore we now need to accomodate the symmetric system
			MatrixOpt::Base_MatrixOpt matrix_opt;
			value = matrix_opt.Sparse_matrix_dot_Vector(this->S_half,value);
	}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::s_r(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::s_phi(const double r,const double phi)const
		{

			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::thetaP(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::R_rr(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::R_thetatheta(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::R_rtheta(const double r,const double phi)const
		{
			return 0;

		}

	template<int dim>
		double
		PoissonHeat_HigherOrder<dim>::R_zz(const double r,const double phi)const
		{
			return 0;

		}

	}