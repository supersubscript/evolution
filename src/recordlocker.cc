#include <iostream>
#include <sstream>
#include <limits>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <unistd.h>
#include <fcntl.h>
#include <map>
#include <ext/stdio_filebuf.h>

#include "recordlocker.h"


recordlocker::recordlocker(const string &dir):
	_dir(dir), _rec(0), _recfd(-1),
	_maxrec(std::numeric_limits<record_t>::max())
{
}

recordlocker::recordlocker(const string &dir, record_t max):
	_dir(dir), _rec(0), _recfd(-1), _maxrec(max)
{
}

recordlocker::~recordlocker()
{
	if(_recfd >= 0)
		close(_recfd);
	// Leave the lock file as a sign that we didn't exit cleanly.
}


// Returns the number of the currently locked record
recordlocker::record_t recordlocker::get_record() const
{
	if(_recfd < 0)
		throw string("state error in get_record: no record is locked");
	return _rec;
}

// Returns the value associated with the currently locked record
recordlocker::value_t recordlocker::get_value() const
{
	if(_recfd < 0)
		throw string("state error in get_value: no record is locked");
	return _value;
}

recordlocker::string recordlocker::recordfilename(__attribute__((unused)) record_t rec) const
{
/*	std::ostringstream s;
	s << _dir << "/record_" << rec << ".lock";
	return s.str();*/
	return _dir + "/records.lock";
}

recordlocker::string recordlocker::recordsfilename() const
{
	return _dir + "/records";
}

static int lock_records(const std::string &fname)
{
	int fd = open(fname.c_str(), O_RDWR | O_CREAT, 0666);
	if(fd < 0)
		throw "failed to open " + fname;

	struct flock fl;
	fl.l_type = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start = fl.l_len = 0;

	if(fcntl(fd, F_SETLKW, &fl) == -1)
	{
		close(fd);
		throw "failed to lock " + fname;
	}
	return fd;
}

static int try_lock_record(const std::string &fname, recordlocker::record_t rec)
{
	int fd = open(fname.c_str(), O_RDWR | O_CREAT, 0666);
	if(fd < 0)
		throw "failed to open " + fname + ": " + strerror(errno);

	struct flock fl;
	fl.l_type = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start = rec;
	fl.l_len = 1;

	if(fcntl(fd, F_SETLK, &fl))
	{
		if(errno == EACCES || errno == EAGAIN)
		{
			close(fd);
			return -1;
		}
		int e = errno;
		close(fd);
		throw "failed to lock " + fname + ": " + strerror(e);
	}
	return fd;
}

void recordlocker::find_free_record(value_t endvalue, bool done)
{
	// We must read all the records before finding a free one, because
	// they may be in random order in the file.
	int recsfd = lock_records(recordsfilename());
	__gnu_cxx::stdio_filebuf<char> fdbuf(recsfd, std::ios::in);
	std::istream in(&fdbuf);
	if(!in)
	{
		close(recsfd);
		throw string("failed to open records as istream");
	}
	std::map<record_t, value_t> records;
	string s;
	while(getline(in, s))
	{
		std::istringstream ss(s);
		record_t rec;
		value_t val;
		ss >> rec >> val;
		if(!ss || !ss.eof())
			throw "record format error on '" + s + "'in " + recordsfilename();
		records[rec] = val;
	}
	if(!in.eof())
		throw "error reading records from file " + recordsfilename();

//	_rec = 0;
	_value = 0;
	while(_rec < _maxrec)
	{
		++_rec;
		auto iter = records.find(_rec);
		if(iter == records.end())
		{
			if(done)
				break;
			_value = 0;
		}
		else
			_value = iter->second;
		if((_value < endvalue) == done)
			continue;

		_recfd = try_lock_record(recordfilename(_rec), _rec);
		if(_recfd >= 0)
		{
			fdbuf.close();
			return;
		}
	}
	// We've done the last record already.
	_rec = 0;

	// This closes recsfd.
	fdbuf.close();
}


recordlocker::record_t recordlocker::lock_record(value_t endvalue)
{
	if(_recfd >= 0)
		throw string("state error in get_record: a record is already locked");
	find_free_record(endvalue, false);
	return _rec;
}

recordlocker::record_t recordlocker::lock_finished_record(value_t endvalue)
{
	if(_recfd >= 0)
		throw string("state error in get_f_record: a record is already locked");
	find_free_record(endvalue, true);
	return _rec;
}

void recordlocker::unlock_record(value_t value)
{
	if(_recfd < 0)
		throw string("state error in unlock_record: no record is locked");
	if(value != _value)
	{
		std::ostringstream ss;
		ss << _rec << "\t" << value << "\n";
		string s = ss.str();
		int recsfd = lock_records(recordsfilename());
		close(_recfd);
		lseek(recsfd, 0, SEEK_END);
		if(write(recsfd, s.c_str(), s.size()) != (int)s.size())
		{
			close(recsfd);
			throw string("failed to write record number to records file");
		}
		close(recsfd);
	}
	else
	{
		// Unlink before closing to avoid race condition if someone picks
		// up this record between the two operations.
		//	unlink(recordfilename(_rec).c_str());
		close(_recfd);
	}

	_recfd = -1;
}
