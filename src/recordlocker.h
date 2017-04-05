#ifndef FOOSIM_RECORDLOCKER_H
#define FOOSIM_RECORDLOCKER_H

#include <string>
#include <limits>

// A class for obtaining work units to be processed, via lock files on
// a shared file system (e.g. NFS). The work units are called 'records'
// here, and each record has a number and an associated value which may
// be used to e.g. store the iteration count in an optimization algortithm.
class recordlocker
{
	typedef std::string string;

  public:
  	// The type of the record number
	typedef uint32_t record_t;
	// The type of the associated value
	typedef uint64_t value_t;

	// Set up a recordlocker with no limit to the number of records
	recordlocker(const string &dir);
	// Set up a recordlocker with a limited number of records
	recordlocker(const string &dir, record_t max);

	// No copy constructor because we keep a lock.
	recordlocker(const recordlocker &) = delete;
	recordlocker & operator=(const recordlocker &) = delete;

	~recordlocker();

	// Returns the number of an available record whose number is less than or
	// equal to the maximum set in the constructor and whose value is less
	// than endvalue. If no such record exists, 0 is returned. On success,
	// unlock_record must then be called before a new record can be locked.
	record_t lock_record(value_t endvalue =
		std::numeric_limits<value_t>::max());

	// Returns the number of an available record whose number is less than or
	// equal to the maximum set in the constructor and whose value is greater
	// than or equal to endvalue. If no such record exists, 0 is returned. On
	// success, unlock_record must then be called before a new record can be
	// locked.
	record_t lock_finished_record(value_t endvalue);

	// Returns the number of the currently locked record
	record_t get_record() const;

	// Returns the value associated with the currently locked record
	value_t get_value() const;

	// Closes the current record, setting its value as per the argument. Thus
	// if you set a higher value, the record can be considered by lock_record
	// as completed.
	void unlock_record(value_t value);

	// Closes the current record without altering its value.
	inline void unlock_record()
	{
		unlock_record(_value);
	}


  private:
	string recordfilename(record_t rec) const;
	string recordsfilename() const;
	void find_free_record(value_t endvalue, bool done);

	string _dir;
	record_t _rec;
	value_t _value;
	int _recfd;
	record_t _maxrec;
};

#endif
