using ErrorAnalysis
using Base.Test

# write your own tests here
@testset "Generic" begin
    x = 3.1234;
    y = 3.021;
    @test iswithinerrorbars(x,y,0.1,true) == false
    @test iswithinerrorbars(x,y,0.11,true) == true
    A = rand(2,2);
    B = rand(2,2);
end
