using ErrorAnalysis
using Base.Test

# write your own tests here
@testset "Generic" begin

    x = 3.1234;
    y = 3.021;
    A = [0.531844 0.717453; 0.552965 0.421109];
    B = [0.273785 0.212329; 0.175248 0.598923];
    @test !iswithinerrorbars(x,y,0.1,true)
    @test iswithinerrorbars(x,y,0.11,true) == true
    @test !iswithinerrorbars(A,B,fill(0.01,size(x)...), true)
    @test iswithinerrorbars(A,B,fill(0.6,size(x)...), true)

end
