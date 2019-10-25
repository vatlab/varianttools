// Copyright 2010 Complete Genomics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License. You
// may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/Exception.hpp"

#include <cerrno>
#include <cstdio>                                 // SEEK_SET, etc.
#include <boost/config.hpp>                       // BOOST_JOIN
#include <boost/iostreams/detail/error.hpp>
#include <boost/iostreams/detail/config/dyn_link.hpp>
#include <boost/iostreams/detail/config/rtl.hpp>  // BOOST_IOSTREAMS_FD_*
#include <boost/iostreams/detail/config/windows_posix.hpp>
#include <boost/iostreams/detail/ios.hpp>         // openmodes, failure.
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/integer_traits.hpp>
#include <boost/system/error_code.hpp>

    // OS-specific headers for low-level i/o.

#include <fcntl.h>       // file opening flags.
#include <sys/stat.h>    // file access permissions.
#ifdef CGA_USE_WIN_API
# include <io.h>         // low-level file i/o.
# define WINDOWS_LEAN_AND_MEAN
# include <windows.h>
# ifndef INVALID_SET_FILE_POINTER
#  define INVALID_SET_FILE_POINTER ((DWORD)-1)
# endif
#else
# include <sys/types.h>  // mode_t.
# include <unistd.h>
# include <string.h>
# include <errno.h>
#endif

namespace cgatools { namespace util {

    //-----------------------------------------------------------------------
    // Devices
    //-----------------------------------------------------------------------

    std::string formatErrorMessage(const std::string& fn, const char* msg)
    {
        using boost::system::error_code;
        // changed by Bo Peng from get_system_category() to system_category() to support newer version of boost
        // using boost::system::system_category;

        std::string buf = std::string(msg) + " for file '" + fn; //  + "': ";
    //#if defined(CGA_USE_WIN_API)
        // changed by Bo Peng from get_system_category() to system_category() to support newer version of boost
        //buf += error_code( ::GetLastError(), system_category() ).message();
    //#else
        // changed by Bo Peng from leaf() to filename() to support newer version of boost
        //buf += error_code( errno, system_category() ).message();
    //#endif
        return buf;
    }

    FileDescriptorDevice::FileDescriptorDevice() : pimpl_(new impl) { }

    FileDescriptorDevice::FileDescriptorDevice(handle_type fd, bool close_on_exit)
        : pimpl_(new impl(fd, close_on_exit))
        { }

    FileDescriptorDevice::FileDescriptorDevice( const char* path,
                                    BOOST_IOS::openmode mode,
                                    BOOST_IOS::openmode base_mode )
        : pimpl_(new impl)
    { open(std::string(path), mode, base_mode); }

    FileDescriptorDevice::FileDescriptorDevice( const std::string& path,
                                    BOOST_IOS::openmode mode,
                                    BOOST_IOS::openmode base_mode )
        : pimpl_(new impl)
    { open(path, mode, base_mode); }


    void FileDescriptorDevice::open( const std::string& path,
                                     BOOST_IOS::openmode m,
                                     BOOST_IOS::openmode base )
    {
        using namespace std;

        pimpl_->fn_ = path;

        m |= base;
        m &= ~BOOST_IOS::binary; // binary is implied for this class

        // ios::ate not supported
        const BOOST_IOS::openmode ALLOWED_FLAGS =
            BOOST_IOS::in | BOOST_IOS::out |
            BOOST_IOS::trunc | BOOST_IOS::app;

        CGA_ASSERT( 0 == (m & ~ALLOWED_FLAGS) );

        bool forRead = (m & BOOST_IOS::in) != 0;
        bool forWrite = (m & BOOST_IOS::out) != 0;
        bool forReadWrite = forRead && forWrite;
        bool truncate = (m & BOOST_IOS::trunc) != 0;
        bool append = (m & BOOST_IOS::app) != 0;

        // at least one of in/out must be specified
        CGA_ASSERT(forRead || forWrite);

        // append can only be combined with out
        CGA_ASSERT(!append || (forWrite && !forRead && !truncate));

        // truncate is only allowed for in+out or out modes
        CGA_ASSERT(!truncate || forWrite);

        // plain ios::out mode implies ios::trunc
        if (forWrite && !forRead && !append) {
            truncate = true;
        }

        if (append) {
            pimpl_->flags_ |= impl::append;
        }

        if ( forWrite ) {
            pimpl_->flags_ |= impl::flushable;
        }

    #ifdef CGA_USE_WIN_API
        DWORD dwDesiredAccess;
        DWORD dwCreationDisposition;
        if (forReadWrite) {
            dwDesiredAccess = GENERIC_READ | GENERIC_WRITE;
            dwCreationDisposition = OPEN_EXISTING;
        } else if (forRead) {
            dwDesiredAccess = GENERIC_READ;
            dwCreationDisposition = OPEN_EXISTING;
        } else if (forWrite) {
            dwDesiredAccess = GENERIC_WRITE;
            dwCreationDisposition = OPEN_ALWAYS;
            // windows doesn't have native append mode; it will be simulated
            // by seeking for every write
        }
        if (truncate) {
            dwCreationDisposition = CREATE_ALWAYS;
        }

        HANDLE handle =
            ::CreateFileA( path.c_str(),
                        dwDesiredAccess,
                        FILE_SHARE_READ | FILE_SHARE_WRITE,
                        NULL,                   // lpSecurityAttributes
                        dwCreationDisposition,
                        FILE_ATTRIBUTE_NORMAL,
                        NULL );                 // hTemplateFile
        if (handle != INVALID_HANDLE_VALUE) {
            pimpl_->handle_ = handle;
            pimpl_->flags_ |= impl::close_on_exit;
        } else {
            pimpl_->flags_ = 0;
            throw Exception(formatErrorMessage(pimpl_->fn_, "open failed"));
        }
    #else

        int oflag = 0;
        if (forReadWrite) {
            oflag |= O_RDWR;
        } else if (forRead) {
            oflag |= O_RDONLY;
        } else if (forWrite) {
            oflag |= O_WRONLY;
            if (append)
                oflag |= O_APPEND;
        }
        if (truncate) {
            oflag |= (O_CREAT|O_TRUNC);
        }

    #ifdef _LARGEFILE64_SOURCE
        oflag |= O_LARGEFILE;
    #endif

        mode_t pmode = S_IRUSR | S_IWUSR |
                    S_IRGRP | S_IWGRP |
                    S_IROTH | S_IWOTH;

        int fd = BOOST_IOSTREAMS_FD_OPEN(path.c_str(), oflag, pmode);
        if (fd == -1) {
            pimpl_->flags_ = 0;
            throw Exception(formatErrorMessage(pimpl_->fn_, "open failed"));
        } else {
            pimpl_->handle_ = fd;
            pimpl_->flags_ |= impl::close_on_exit;
        }
    #endif


    }

    void FileDescriptorDevice::open
        ( const char* path, BOOST_IOS::openmode m,
        BOOST_IOS::openmode base )
    { open(std::string(path), m, base); }

    std::streamsize FileDescriptorDevice::read(char_type* s, std::streamsize n)
    {
    #ifdef CGA_USE_WIN_API
        DWORD result;
        if (!::ReadFile(pimpl_->handle_, s, n, &result, NULL))
            throw Exception(formatErrorMessage(pimpl_->fn_, "read failed"));
        return result == 0 ? -1 : static_cast<std::streamsize>(result);
    #else
        errno = 0;
        std::streamsize result = BOOST_IOSTREAMS_FD_READ(pimpl_->handle_, s, n);
        if (errno != 0)
            throw Exception(formatErrorMessage(pimpl_->fn_, "read failed"));
        return result == 0 ? -1 : result;
    #endif
    }

    std::streamsize FileDescriptorDevice::write(const char_type* s, std::streamsize n)
    {
    #ifdef CGA_USE_WIN_API
        if (pimpl_->flags_ & impl::append) {
            DWORD const dwResult =
                ::SetFilePointer(pimpl_->handle_, 0, NULL, FILE_END);
            if ( dwResult == INVALID_SET_FILE_POINTER &&
                ::GetLastError() != NO_ERROR )
            {
                throw Exception(formatErrorMessage(pimpl_->fn_, "append seek failed"));
            }
        }
        DWORD ignore;
        if (!::WriteFile(pimpl_->handle_, s, n, &ignore, NULL))
            throw Exception(formatErrorMessage(pimpl_->fn_, "write failed"));
        return n;
    #else
        int amt = BOOST_IOSTREAMS_FD_WRITE(pimpl_->handle_, s, n);
        if (amt < n)
            throw Exception(formatErrorMessage(pimpl_->fn_, "write failed"));
        return n;
    #endif
    }

    std::streampos FileDescriptorDevice::seek(boost::iostreams::stream_offset off,
                                        BOOST_IOS::seekdir way)
    {
        using namespace std;
    #ifdef CGA_USE_WIN_API
        LONG lDistanceToMove = static_cast<LONG>(off & 0xffffffff);
        LONG lDistanceToMoveHigh = static_cast<LONG>(off >> 32);
        DWORD dwResultLow =
            ::SetFilePointer( pimpl_->handle_,
                            lDistanceToMove,
                            &lDistanceToMoveHigh,
                            way == BOOST_IOS::beg ?
                                FILE_BEGIN :
                                way == BOOST_IOS::cur ?
                                    FILE_CURRENT :
                                    FILE_END );
        if ( dwResultLow == INVALID_SET_FILE_POINTER &&
            ::GetLastError() != NO_ERROR )
        {
            throw Exception(formatErrorMessage(pimpl_->fn_, "seek failed"));
        } else {
            return boost::iostreams::offset_to_position(
                    (boost::iostreams::stream_offset(lDistanceToMoveHigh) << 32) + dwResultLow
                );
        }
    #else
        if ( off > boost::integer_traits<BOOST_IOSTREAMS_FD_OFFSET>::const_max ||
            off < boost::integer_traits<BOOST_IOSTREAMS_FD_OFFSET>::const_min )
        {
            throw BOOST_IOSTREAMS_FAILURE("bad offset");
        }
        boost::iostreams::stream_offset result =
            BOOST_IOSTREAMS_FD_SEEK(
                pimpl_->handle_,
                static_cast<BOOST_IOSTREAMS_FD_OFFSET>(off),
                ( way == BOOST_IOS::beg ?
                    SEEK_SET :
                    way == BOOST_IOS::cur ?
                        SEEK_CUR :
                        SEEK_END )
            );
        if (result == -1)
            throw Exception(formatErrorMessage(pimpl_->fn_, "seek failed"));
        return boost::iostreams::offset_to_position(result);
    #endif
    }

    void FileDescriptorDevice::close()
    {
        close_impl(*pimpl_);
    }

    bool FileDescriptorDevice::fsync()
    {
        fsync_impl(*pimpl_);
        return true;
    }

    void FileDescriptorDevice::fsync_impl(impl& i)
    {
        if (i.flags_ & impl::flushable) {
    #ifdef CGA_USE_WIN_API
            if (! ::FlushFileBuffers(i.handle_) ) {
                throw Exception(formatErrorMessage(i.fn_, "fsync failed"));
            }
    #else
            int rc = ::fsync(i.handle_);
            if (rc == -1) {
                int ec = errno;
                if (ec == EROFS || ec == EINVAL) {
                    // According to fsync(2) man page, these error codes indicate
                    // that the descriptor is bound to a special device that does
                    // not support sync. In practice it happens if the file is
                    // actually /dev/null. Silently ignore the error.
                } else {
                    throw Exception(formatErrorMessage(i.fn_, "fsync failed"));
                }
            }
    #endif
        }
    }

    void FileDescriptorDevice::close_impl(impl& i)
    {
    #ifdef CGA_USE_WIN_API
        if (i.handle_ != reinterpret_cast<handle_type>(-1)) {
            try {
                fsync_impl(i);
            } catch (const std::exception&) {
                ::CloseHandle(i.handle_);
                i.handle_ = reinterpret_cast<handle_type>(-1);
                i.flags_ = 0;
                throw;
            }

            BOOL result = ::CloseHandle(i.handle_);
            i.handle_ = reinterpret_cast<handle_type>(-1);
            i.flags_ = 0;

            if (!result) {
                throw Exception(formatErrorMessage(i.fn_, "close failed"));
            }
        }
    #else
        if (i.handle_ != -1) {
            try {
                fsync_impl(i);
            } catch (const std::exception&) {
                BOOST_IOSTREAMS_FD_CLOSE(i.handle_);
                i.handle_ = -1;
                i.flags_ = 0;
                throw;
            }

            int result = BOOST_IOSTREAMS_FD_CLOSE(i.handle_);
            i.handle_ = -1;
            i.flags_ = 0;

            if (result == -1) {
                throw Exception(formatErrorMessage(i.fn_, "close failed"));
            }
        }
    #endif
    }

    //-----------------------------------------------------------------------
    // Input/OutputStream
    //-----------------------------------------------------------------------

    void InputStream::open(const char* fn)
    {
        base_type::exceptions(std::ios::badbit);
        base_type::open(FileSourceDevice(fn));
    }

    bool InputStream::isReadable(const std::string& fn)
    {
        try
        {
            InputStream in(fn);
        }
        catch(const std::exception&)
        {
            return false;
        }
        return true;
    }

    boost::shared_ptr<std::istream>
    InputStream::openCompressedInputStreamByExtension(const std::string& fn)
    {
        boost::shared_ptr<std::istream> result;
        if (fn.size() >= 3 && fn.compare(fn.size() - 3, 3, ".gz") == 0)
            result.reset(new CompressedInputStream(fn));
        else if (fn.size() >= 4 && fn.compare(fn.size() - 4, 4, ".bz2") == 0)
        {
            CompressedInputStream* in = new CompressedInputStream();
            result.reset(in);
            in->openBZ2(fn);
        }
        else
            result.reset(new InputStream(fn));
        return result;
    }

    bool InputStream::getline(std::istream& in, std::string& line)
    {
        std::getline(in, line);
        if (line.size() > 0 && '\r' == line[line.size()-1])
            line.resize(line.size()-1);
        return in.good();
    }

    void OutputStream::open(const char* fn)
    {
        base_type::exceptions(std::ios::badbit | std::ios::failbit);
        base_type::open(FileSinkDevice(fn));
    }

    OutputStream::~OutputStream()
    {
        try {
            // explicitly call close() since Boost's destructor would
            // catch and ignore exceptions from it, and we want to
            // kill the application instead
            //if (is_open()) {
                close(); // close() checks is_open() itself
            //}
        } catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            CGA_ASSERT(false);
        }
    }

    //-----------------------------------------------------------------------
    // CompressedInput/OutputStream
    //-----------------------------------------------------------------------
    void CompressedInputStream::open(const char* fn)
    {
        base_type::push(boost::iostreams::gzip_decompressor(), 4*1024);
        base_type::push(FileSourceDevice(fn));
        base_type::exceptions(std::ios::badbit);
    }

    void CompressedInputStream::openBZ2(const char* fn)
    {
        base_type::push(boost::iostreams::bzip2_decompressor(), 4*1024);
        base_type::push(FileSourceDevice(fn));
        base_type::exceptions(std::ios::badbit);
    }

    void CompressedOutputStream::open(const char* fn, int clev)
    {
        base_type::push(boost::iostreams::gzip_compressor(clev), 4*1024);
        base_type::push(FileSinkDevice(fn));
        base_type::exceptions(std::ios::badbit | std::ios::failbit);
    }

    void CompressedOutputStream::openBZ2(const char* fn)
    {
        base_type::push(boost::iostreams::bzip2_compressor(), 4*1024);
        base_type::push(FileSinkDevice(fn));
        base_type::exceptions(std::ios::badbit | std::ios::failbit);
    }

    boost::shared_ptr<std::ostream>
    OutputStream::openCompressedOutputStreamByExtension(const std::string& fn)
    {
        boost::shared_ptr<std::ostream> result;
        if (fn.size() >= 3 && fn.compare(fn.size() - 3, 3, ".gz") == 0)
            result.reset(new CompressedOutputStream(fn));
        else if (fn.size() >= 4 && fn.compare(fn.size() - 4, 4, ".bz2") == 0)
        {
            CompressedOutputStream* out = new CompressedOutputStream();
            result.reset(out);
            out->openBZ2(fn);
        }
        else
            result.reset(new OutputStream(fn));
        return result;
    }

    CompressedOutputStream::~CompressedOutputStream()
    {
        try {
            // pop and close all elements of the pipeline
            reset();
        } catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            CGA_ASSERT(false);
        }
    }

    void writeBinaryBool(std::ostream& out, bool val)
    {
        char cVal = val ? 1 : 0;
        out.put(cVal);
    }

    void writeBinaryString(std::ostream& out, const std::string& val)
    {
        writeBinaryUIntZC(out, val.size());
        out.write(val.c_str(), val.size());
    }

    void readBinaryBool(std::istream& in, bool* val)
    {
        char ch;
        in.get(ch);
        if (!in.good())
            throw Exception("failed to read binary bool: unexpected eof");
        *val = 0 != ch;
    }

    void readBinaryString(std::istream& in, std::string* val)
    {
        size_t sz;
        readBinaryUIntZC(in, &sz);
        val->resize(sz);
        in.read(&(*val)[0], val->size());
        if (!in.good())
            throw Exception("failed to read binary string: unexpected eof");
    }

} } // cgatools::util
