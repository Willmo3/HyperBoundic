//
// Created by will on 9/11/25.
//

#ifndef PDEAPPROX_SYSTEMAPPROXIMATION_H
#define PDEAPPROX_SYSTEMAPPROXIMATION_H
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>

#include "solvers/flux/FluxFunction.hpp"
#include "domains/Numeric.hpp"
#include "CflCheck.hpp"

#include "cereal/archives/json.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/unordered_map.hpp"
#include "Waffine/WaffineForm.hpp"

/**
 * Discretization of a physical system represented by a hyperbolic PDE.
 * This mesh uses a fixed spatial dimension, ideal for finite difference methods.
 * @param T numeric type to approximate system.
 */
template<typename T>
requires Numeric<T>
class RectangularMesh {
public:
    /*
     * Constructors
     */

    /**
     * Empty discretization matrix.
     * @param discretization_size Number of spatial discretization points, > 0.
     * @param num_timesteps Number of timesteps for this discretization->
     */
    RectangularMesh(uint32_t discretization_size, uint32_t num_timesteps)
        :_discretization_size(discretization_size),  _num_timesteps(num_timesteps) {
        assert(discretization_size > 0);
        assert(num_timesteps > 0);

        _system = static_cast<T *>(calloc(sizeof(T), num_timesteps * discretization_size));
        assert(_system);
    }
    /**
     * Copy initial conditions into discretization matrix.
     * @param initial_conditions Array of starting conditions for the system, of len discretization_size.
     * We manually assign each element to ensure copy constructors are called.
     */
    void copy_initial_conditions(const std::vector<T> &initial_conditions) {
        assert(initial_conditions.size() == discretization_size());
        for (auto timestep = 0; timestep < _num_timesteps; timestep++) {
            _system[timestep] = initial_conditions[timestep];
        }
    }
    /**
     * Destructor
     */
    ~RectangularMesh() {
        free(_system);
        _system = nullptr;
    }

    /*
     * Accessors
     */
    uint32_t discretization_size() const {
        return _discretization_size;
    }
    uint32_t num_timesteps() const {
        return _num_timesteps;
    }
    T get(uint32_t timestep, uint32_t index) const {
        assert(timestep < _num_timesteps);
        assert(index < _discretization_size);
        return _system[timestep * _discretization_size + index];
    }
    void set(uint32_t timestep, uint32_t index, T value) {
        assert(timestep < _num_timesteps);
        assert(index < _discretization_size);
        _system[timestep * _discretization_size + index] = value;
    }

    /*
     * Serialization
     */

    /**
     * @return A json representation of this data.
     */
    std::string to_json_string() {
        std::ostringstream ss;
        // Inner scope needed to ensure proper flushing.
        {
            cereal::JSONOutputArchive o(ss);
            o(FixedMeshData(this));
            ss.flush();
        }
        return ss.str();
    }

    /**
     * @param strrep JSON string representation of a discretization.
     * @return A new discretization object created from this string.
     */
    static RectangularMesh from_json_string(const std::string &strrep) {
        auto input = std::istringstream(strrep);
        FixedMeshData data_tuple = FixedMeshData();
        // Place in inner scope to ensure proper flushing.
        {
            cereal::JSONInputArchive archive(input);
            archive(data_tuple);
        }

        // Manual deep copy to handle unordered maps in affine forms.
        // NOTE: simple memcpy does not cut it -- we need to call copy constructor!
        auto data = static_cast<T *>(calloc(data_tuple.discretization_size * data_tuple.num_timesteps, sizeof(T)));
        for (auto i = 0; i < data_tuple.discretization_size * data_tuple.num_timesteps; i++) {
            data[i] = data_tuple.system[i];
        }

        return RectangularMesh(data_tuple.discretization_size, data_tuple.num_timesteps, data);
    }

    /*
     * Assorted helpers
     */

    void print_system() const {
        for (auto t = 0; t < _num_timesteps; t++) {
            std::cout << "T" << t << ": ";
            for (auto i = 0; i < _discretization_size; i++) {
                std::cout << _system[_discretization_size * t + i] << " ";
            }
            std::cout << std::endl;
        }
    }

    bool equals(const RectangularMesh &other) const {
        if (other._discretization_size != _discretization_size || other._num_timesteps != _num_timesteps) {
            return false;
        }

        for (auto i = 0; i < _discretization_size * _num_timesteps; i++) {
            if (other._system[i] != _system[i]) {
                return false;
            }
        }

        return true;
    }

private:
    // Using raw pointer to enable low-level mem management -- i.e. transfer to GPU
    T *_system;
    const uint32_t _discretization_size;
    const uint32_t _num_timesteps;

    /**
     * Internal data tuple for reading/writing
     */
    struct FixedMeshData {
        uint32_t discretization_size;
        uint32_t num_timesteps;
        std::vector<T> system;

        /**
         * @param discretization Discretization to construct into serializable form.
         */
        explicit FixedMeshData(const RectangularMesh *discretization):
            discretization_size(discretization->discretization_size()), num_timesteps(discretization->num_timesteps()),
            system(std::vector<T>(discretization->discretization_size() * discretization->num_timesteps())) {

            for (auto i = 0; i < discretization->discretization_size() * discretization->num_timesteps(); i++) {
                system[i] = discretization->_system[i];
            }
        }

        /**
         * Dummy constructor for reading in deserialized values.
         */
        FixedMeshData(): discretization_size(0), num_timesteps(0), system(std::vector<T>()) {}

        /*
         * Serialization support
         */
        template<class Archive>
        void serialize(Archive & archive) {
            archive(cereal::make_nvp("discretization_size", discretization_size),
                cereal::make_nvp("num_timesteps", num_timesteps),
                cereal::make_nvp("system", system));
        }
    };

    /**
     * Explicit value constructor -- internal only!
     *
     * @param discretization_size Size of discretization, > 0.
     * @param num_timesteps Number of timesteps for this discretization->
     * @param system Array of starting conditions for the system, of len discretization_size.
     */
    RectangularMesh(uint32_t discretization_size, uint32_t num_timesteps, T *system):
        _discretization_size(discretization_size), _num_timesteps(num_timesteps), _system(system) {
        assert(system);
    }
};

#endif //PDEAPPROX_SYSTEMAPPROXIMATION_H