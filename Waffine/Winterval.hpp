//
// Created by will on 9/15/25.
//

#ifndef WINTERVAL_WINTERVAL_H
#define WINTERVAL_WINTERVAL_H
#include <string>

/**
 * @brief Intervals defined with support for operator overloading.
 * Intended for use in simple numeric abstract interpretation.
 *
 * @author Will Morris
 */
class Winterval {
public:
    /**
     * @param min minimum value in interval, inclusive.
     * @param max Maximum value in interval, inclusive.
     */
    Winterval(double min, double max);
    ~Winterval();

    /*
     * Accessors
     */
    double min() const;
    double max() const;
    /**
     * @return the midpoint of this interval.
     */
    double mid() const;
    /**
     * @return the distance of the min and max to the midpoint.
     */
    double radius() const;

    /*
     * Binary Winterval operations.
     * Note: in order to enable nested operator application, we don't use references.
     * Overhead should be small, but if we need to squeeze out some cycles later, look here!
     */
    Winterval operator+(const Winterval &rhs) const;
    Winterval operator-(const Winterval &rhs) const;
    Winterval operator*(const Winterval &rhs) const;
    Winterval operator/(const Winterval &rhs) const;

    /*
     * Relational Winterval operations.
     * Useful for bounds checks in abstract interpretation.
     */
    bool operator==(const Winterval &rhs) const;
    /*
     * Threshold operators.
     */
    // Evaluated against max bound.
    bool operator<=(double rhs) const;
    bool operator<(double rhs) const;
    // Evaluated against min bound.
    bool operator>=(double rhs) const;
    bool operator>(double rhs) const;

    /*
     * Scalar Winterval operations.
     */
    Winterval operator+(double rhs) const;
    Winterval operator-(double rhs) const;
    Winterval operator*(double rhs) const;
    Winterval operator/(double rhs) const;

    /*
     * Other numeric functions.
     */
    // TODO: as more flux functions are needed, add them here.
    /**
     * Compute the hyperbolic tangent of an interval.
     * Since tanh is monotonic, it is sufficient to apply it to the lower, upper bounds (i.e. optimal conditions on boundary)
     * @return A new interval whose value is the hyperbolic tangent of this interval.
     */
    Winterval tanh() const;

    /**
     * Compute an interval to some constant power.
     * @param power Power to raise the interval to.
     * Note that for intervals, integer powers should always be nonnegative -- we represent as an int to comply with idiomatic Cpp.
     *
     * @return a new interval equal to this one raised to power
     */
    Winterval pow(int power) const;

    /**
     * @return The absolute value of this interval.
     */
    Winterval abs() const;

    /**
     * @brief Inclusive bounds check.
     * @param value Value to check if in interval
     * @return Whether value is in the interval.
     */
    bool contains(double value) const;
private:
    double _min;
    double _max;
};

std::ostream& operator<<(std::ostream& os, Winterval rhs);

#endif //WINTERVAL_WINTERVAL_H