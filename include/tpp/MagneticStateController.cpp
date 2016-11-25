// Constructor
template<class INTE>
MagneticStateController<INTE>::MagneticStateController(
    const ParticleCluster _geom,
    const maggie:: field _happ,
    const StocLLG::vector _llgs,
    const std::vector<maggie::moment> init_state,
    const double time_step,
    const unsigned int seed
    )
    : geom( _geom )
    , happ( _happ )
    , llgs( _llgs )
    , N( geom.getNParticles() )
    , disp( geom.getReducedDistancesRef() )
    , unit_disp( boost::extents[N][N][3] )
    , cube_disps( boost::extents[N][N] )
    , v( geom.getReducedVolumes() )
    , k( geom.getReducedAnisConstants() )
{
    // Set RNG seed
    integrator_rng.seed( seed );

    // Compute dipole-dipole interaction prefactor
    ddNorm = Constants::MU0 * geom.getParticle(0).getMs()
        * geom.getParticle(0).getMs()
        / ( 8 * M_PI * geom.getAverageAnisConstant() );

    // Compute the cube distances and displacement unit vectors
    for( int i=0; i!=N; ++i )
        for( int j=0; j!=N; ++j )
        {
            double mag = sqrt(
                disp[i][j][0] * disp[i][j][0]
                + disp[i][j][1] * disp[i][j][1]
                + disp[i][j][2] * disp[i][j][2]
                );
            unit_disp[i][j][0] = disp[i][j][0] / mag;
            unit_disp[i][j][1] = disp[i][j][1] / mag;
            unit_disp[i][j][2] = disp[i][j][2] / mag;
            cube_disps[i][j] = mag * mag * mag;
        }

    // Create a vector of field references
    for( int i=0; i!=N; ++i )
        heff.push_back( &(llgs[i].getReducedFieldRef()) );

    // Create the integrators
    for( int i=0; i!=N; ++i )
    {
        integrators.push_back(
            INTE(
                llgs[i], init_state[i], 0.0, time_step, integrator_rng
                )
            );
    }
    for( int i=0; i!=N; ++i )
        states.push_back( &( integrators[i].getState() ) );
}

template<class INTE>
void MagneticStateController<INTE>::updateEffectiveField()
{
    // For each particle in the ensemble
    for( int i=0; i<N; ++i )
    {
        // Get the current state
        moment state = integrators[i].getState();

        // Get the effective field reference
        field& heff_ref = *(heff[i]);

        // Initialise the state with the anisotropy contribution
        geom.getParticle( i ).computeAnisotropyField(heff_ref, k[i], state);


        // Add the Zeeman field contribution
        heff_ref[0] = heff_ref[0] + happ[0];
        heff_ref[1] = heff_ref[1] + happ[1];
        heff_ref[2] = heff_ref[2] + happ[2];

        // Add the contribution from the interacting field
        for( int j=0; j<i; ++j )
        {
            moment coupled_state = integrators[j].getState();
            double dot_prod = coupled_state[0] * unit_disp[i][j][0]
                + coupled_state[1] * unit_disp[i][j][1]
                + coupled_state[2] * unit_disp[i][j][2];

            heff_ref[0] = heff_ref[0] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][0] - coupled_state[0] );
            heff_ref[1] = heff_ref[1] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][1] - coupled_state[1] );
            heff_ref[2] = heff_ref[2] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][2] - coupled_state[2] );
        }
        // Split these loops to avoid continue statement (bad for CPU optimisation)
        for( int j=i+1; j<N; ++j )
        {
            moment coupled_state = integrators[j].getState();
            double dot_prod = coupled_state[0] * unit_disp[i][j][0]
                + coupled_state[1] * unit_disp[i][j][1]
                + coupled_state[2] * unit_disp[i][j][2];

            heff_ref[0] = heff_ref[0] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][0] - coupled_state[0] );
            heff_ref[1] = heff_ref[1] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][1] - coupled_state[1] );
            heff_ref[2] = heff_ref[2] + ddNorm * v[j] / cube_disps[i][j]
                * ( 3 * dot_prod * unit_disp[i][j][2] - coupled_state[2] );
        }
    } // END for each particle in cluster
} // END updateEffectiveField()

template<class INTE>
void MagneticStateController<INTE>::renormaliseStates()
{

    for( int i=0; i<N; ++i )
    {
        double mx=states[i]->at(0);
        double my=states[i]->at(1);
        double mz=states[i]->at(2);

        double norm = sqrt( mx*mx + my*my + mz*mz );
        moment newstate{ mx/norm, my/norm, mz/norm };

        integrators[i].setState( newstate );
    } // END for each particle in the cluster
} // END renormaliseStates()

template<class INTE>
void MagneticStateController<INTE>::step()
{
    updateEffectiveField();
    for( int i=0; i!=N; ++i )
        integrators[i].step();
    renormaliseStates();
}

template<class INTE>
void MagneticStateController<INTE>::reset()
{
    for( unsigned int i=0; i!=N; ++i )
        this->integrators[i].reset();
}

template<class INTE>
void MagneticStateController<INTE>::reset(std::vector<maggie::moment> init_state)
{
    for( unsigned int i=0; i!=N; ++i )
        integrators[i].reset( init_state[i] );
}

template<class INTE>
double MagneticStateController<INTE>::getDDPrefactor() const
{
    return ddNorm;
}

template<class INTE>
matrix_d MagneticStateController<INTE>::getCubeDisplacements() const
{
    return cube_disps;
}

template<class INTE>
const std::vector<moment const*>& MagneticStateController<INTE>::getState() const
{
    return states;
}

template<class INTE>
double MagneticStateController<INTE>::getTime() const
{
    return integrators[0].getTime();
}

template<class INTE>
std::vector<maggie::field> MagneticStateController<INTE>::getEffectiveField() const
{
    std::vector<maggie::field> fields;
    for ( auto heff_ptr : heff )
        fields.push_back( *heff_ptr );
    return fields;
}
