#include <fmt/core.h>
#include <highs/highs.h>    

int main()
{
    fmt::print("Hello World!\n");
    fmt::print("Highs version: {}\n", highsVersion());
    fmt::print("Highs version major: {}\n", highsVersionMajor());
    fmt::print("Highs version minor: {}\n", highsVersionMinor());
    fmt::print("Highs version patch: {}\n", highsVersionPatch());
    fmt::print("Highs githash: {}\n", highsGithash());
    return 0;
}