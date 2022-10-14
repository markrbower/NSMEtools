function [Pxx1,Pxx2,CPxx,MCxx] = Bower_coherence(fname)
    close all
    load( fname );
    figure(1)
    subplot( 2, 1, 1 )
    plot(x1,'r');
    subplot( 2, 1, 2 )
    plot(x2 - 10,'g');
    Fs = 1000;
    Pxx1 = pwelch( x1, [], [], [], Fs );
    Pxx2 = pwelch( x2, [], [], [], Fs );
    CPxx = cpsd( x1, x2, [], [], [], Fs );
    MCxx = mscohere( x1, x2, [], [], [], Fs );
    figure(2)
    subplot(4,1,1)
    plot( log10(Pxx1) );
    subplot(4,1,2)
    plot( log10(Pxx2) );
    subplot(4,1,3)
    cpsd( x1, x2, [], [], [], Fs );
    set( gca, 'YLabel', [] );
    subplot(4,1,4)
    plot( MCxx );
    figure(3)
    [Mxx1,Mxxc1,F] = pmtm( x1, [], [], Fs );
    [Mxx2,Mxxc2,F] = pmtm( x2, [], [], Fs );
    subplot(4,1,1)
    hold on
    plot( log10( Mxx1 ) );
    plot( log10( Mxxc1(:,1) ), 'g' )
    plot( log10( Mxxc1(:,2) ), 'g' )
    subplot(4,1,2)
    hold on
    plot( log10( Mxx1 ) );
    plot( log10( Mxxc1(:,1) ), 'g' )
    plot( log10( Mxxc1(:,2) ), 'g' )
    set( gca, 'YLim', [-2 -.5] );
    set( gca, 'XLim', [0 1E5] );

    subplot(4,1,3)
    hold on
    plot( log10( Mxx2 ) );
    plot( log10( Mxxc2(:,1) ), 'g' )
    plot( log10( Mxxc2(:,2) ), 'g' )
    subplot(4,1,4)
    hold on
    plot( log10( Mxx2 ) );
    plot( log10( Mxxc2(:,1) ), 'g' )
    plot( log10( Mxxc2(:,2) ), 'g' )
    set( gca, 'YLim', [-2 -.5] );
    set( gca, 'XLim', [0 1E5] );

end