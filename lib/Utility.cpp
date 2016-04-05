#include<Utility.hpp>

void progress_bar( std::string taskname, double progress )
{
    int length = 20;
    cout << taskname << "  |";
    int complete = int( 20 * progress );
    for( int i=0; i!=length; ++i )
        i <= complete ? cout << "=" : cout << ".";
    cout << "|\r" << flush;
}
