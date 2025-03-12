#ifndef MATH_UTILS_H
#define MATH_UTILS_H

// Custom clamp function (C++11 compatible for Visual Studio 2013)
template <typename T>
T clamp(T value, T min_value, T max_value)
{
    return (value < min_value) ? min_value : (value > max_value) ? max_value : value;
}

#endif // MATH_UTILS_H

