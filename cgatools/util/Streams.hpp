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

#ifndef CGATOOLS_UTIL_STREAMS_HPP_
#define CGATOOLS_UTIL_STREAMS_HPP_ 1

//! @file Streams.hpp
//! File containing definitions of InputStream, OutputStream,
//! CompressedInputStream, and CompressedOutputStream.

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace cgatools { namespace util {

    // Direct adaptation of boost::iostreams::file_descriptor
    // with better error messages, and an explicit fsync(2) call
    // in close() method
    class  FileDescriptorDevice {
    public:
#ifdef CGA_USE_WIN_API
        typedef void*  handle_type;  // A.k.a HANDLE
#else
        typedef int    handle_type;
#endif
        typedef char   char_type;
        struct category
            : boost::iostreams::seekable_device_tag,
              boost::iostreams::closable_tag
            { };
        FileDescriptorDevice();
        explicit FileDescriptorDevice(handle_type fd, bool close_on_exit = false);

        explicit FileDescriptorDevice( const std::string& path,
                                  BOOST_IOS::openmode mode =
                                      BOOST_IOS::in | BOOST_IOS::out,
                                  BOOST_IOS::openmode base_mode =
                                      BOOST_IOS::in | BOOST_IOS::out );
        explicit FileDescriptorDevice( const char* path,
                                  BOOST_IOS::openmode mode =
                                      BOOST_IOS::in | BOOST_IOS::out,
                                  BOOST_IOS::openmode base_mode =
                                      BOOST_IOS::in | BOOST_IOS::out );
        void open( const std::string& path,
                   BOOST_IOS::openmode =
                       BOOST_IOS::in | BOOST_IOS::out,
                   BOOST_IOS::openmode base_mode =
                       BOOST_IOS::in | BOOST_IOS::out );
        void open( const char* path,
                   BOOST_IOS::openmode =
                       BOOST_IOS::in | BOOST_IOS::out,
                   BOOST_IOS::openmode base_mode =
                       BOOST_IOS::in | BOOST_IOS::out );
        bool is_open() const { return pimpl_->flags_ != 0; }
        std::streamsize read(char_type* s, std::streamsize n);
        std::streamsize write(const char_type* s, std::streamsize n);
        std::streampos seek(boost::iostreams::stream_offset off, BOOST_IOS::seekdir way);
        bool fsync();
        void close();
        handle_type handle() const { return pimpl_->handle_; }
        const std::string& fn() const { return pimpl_->fn_; }
    private:
        struct impl {
            impl() :
#ifdef CGA_USE_WIN_API
                    handle_(reinterpret_cast<handle_type>(-1)),
#else
                    handle_(-1),
#endif
                    flags_(0)
                { }

            impl(handle_type fd, bool close_on_exit)
                : fn_("<handle>"), handle_(fd), flags_(0)
            {
                if (close_on_exit) flags_ |= impl::close_on_exit;
            }

            ~impl()
            {
                if (flags_ & close_on_exit) close_impl(*this);
            }
            enum flags {
                close_on_exit = 1,
                append = 4,
                flushable = 8
            };
            std::string  fn_;
            handle_type  handle_;
            int          flags_;
        };
        friend struct impl;

        static void close_impl(impl&);
        static void fsync_impl(impl&);

        boost::shared_ptr<impl> pimpl_;
    };

    struct FileSourceDevice : private FileDescriptorDevice {
        typedef FileDescriptorDevice::handle_type handle_type;
        typedef char   char_type;
        struct category
          : boost::iostreams::input_seekable,
            boost::iostreams::device_tag,
            boost::iostreams::closable_tag
          { };
        using FileDescriptorDevice::read;
        using FileDescriptorDevice::seek;
        using FileDescriptorDevice::open;
        using FileDescriptorDevice::is_open;
        using FileDescriptorDevice::close;
        using FileDescriptorDevice::handle;
        FileSourceDevice() { }
        explicit FileSourceDevice(handle_type fd, bool close_on_exit = false)
            : FileDescriptorDevice(fd, close_on_exit)
            { }

        explicit FileSourceDevice( const std::string& path,
                                         BOOST_IOS::openmode m = BOOST_IOS::in )
            : FileDescriptorDevice(path, m & ~BOOST_IOS::out, BOOST_IOS::in)
            { }
        explicit FileSourceDevice( const char* path,
                                         BOOST_IOS::openmode m = BOOST_IOS::in )
            : FileDescriptorDevice(path, m & ~BOOST_IOS::out, BOOST_IOS::in)
            { }
    };

    struct FileSinkDevice : private FileDescriptorDevice {
        typedef FileDescriptorDevice::handle_type handle_type;
        typedef char   char_type;
        struct category
          : boost::iostreams::output_seekable,
            boost::iostreams::device_tag,
            boost::iostreams::closable_tag
            // note that the device is NOT marked as flushable;
            // see the comment in OutputStream
          { };
        using FileDescriptorDevice::write;
        using FileDescriptorDevice::seek;
        using FileDescriptorDevice::open;
        using FileDescriptorDevice::is_open;
        using FileDescriptorDevice::close;
        using FileDescriptorDevice::handle;
        using FileDescriptorDevice::fsync;
        FileSinkDevice() { }
        explicit FileSinkDevice(handle_type fd, bool close_on_exit = false)
            : FileDescriptorDevice(fd, close_on_exit)
            { }
        explicit FileSinkDevice( const std::string& path,
                                       BOOST_IOS::openmode m = BOOST_IOS::out )
            : FileDescriptorDevice(path, m & ~BOOST_IOS::in, BOOST_IOS::out)
            { }
        explicit FileSinkDevice( const char* path,
                                       BOOST_IOS::openmode m = BOOST_IOS::out )
            : FileDescriptorDevice(path, m & ~BOOST_IOS::in, BOOST_IOS::out)
            { }
    };

    //! A std::istream with good error messages.
    class InputStream:
        public boost::iostreams::stream<FileSourceDevice>
    {
    public:
        InputStream() {}
        InputStream(const char* fn)
        {
            open(fn);
        }

        InputStream(const std::string& fn)
        {
            open(fn);
        }

        void open(const std::string& fn)
        {
            open(fn.c_str());
        }
        void open(const char* fn);

        void close()
        {
            if (is_open()) {
                base_type::close();
            }
        }

        static bool isReadable(const std::string& fn);

        //! Returns a pointer to an InputStream or
        //! CompressedInputStream, depending on the file extension. File
        //! extensions supported:
        //! - .gz -> gzip format
        //! - .bz2 -> bzip2 format
        //! - other -> straight InputStream
        static boost::shared_ptr<std::istream>
        openCompressedInputStreamByExtension(const std::string& fn);

        //! Retrieves the line from the given input stream, and removes
        //! a carriage return from the end if there is one.
        static bool getline(std::istream& in, std::string& line);

    private:
        typedef boost::iostreams::stream<FileSourceDevice> base_type;
    };

    //! a std::ostream with good error messages, and which does fsync on
    //! close on UNIX.  Note that calling flush() of this class will not cause
    //! fsync(2) call, since FileSinkDevice is not marked as
    //! flushable. This is so that we don't fsync on every f <<
    //! std::endl. If you really need to sync to disk without closing
    //! the stream, you can do:
    //! @code
    //! OutputStream f("foo");
    //! f.flush();
    //! (*f).fsync();
    //! @endcode
    //! This will use operator*() to get access to the underlying
    //! device. This class will fsync automatically when closing
    //! itself.
    class OutputStream:
        public boost::iostreams::stream<FileSinkDevice>
    {
    public:
        OutputStream() {}
        OutputStream(const char* fn)
        {
            open(fn);
        }
        OutputStream(const std::string& fn)
        {
            open(fn);
        }
        ~OutputStream();

        void open(const std::string& fn)
        {
            open(fn.c_str());
        }
        void open(const char* fn);

        void close()
        {
            if (is_open()) {
                base_type::close();
            }
        }

        //! Returns a pointer to an OutputStream or
        //! CompressedOutputStream, depending on the file
        //! extension. File extensions supported:
        //! - .gz -> gzip format
        //! - .bz2 -> bzip2 format
        //! - other -> straight InputStream
        static boost::shared_ptr<std::ostream>
        openCompressedOutputStreamByExtension(const std::string& fn);
    private:
        typedef boost::iostreams::stream<FileSinkDevice> base_type;
    };

    class CompressedInputStream:
        public boost::iostreams::filtering_istream
    {
    public:
        CompressedInputStream() {}
        CompressedInputStream(const char* fn)
        {
            open(fn);
        }
        CompressedInputStream(const std::string& fn)
        {
            open(fn);
        }

        //! Open a gzip-compressed file.
        void open(const std::string& fn)
        {
            open(fn.c_str());
        }
        //! Open a gzip-compressed file.
        void open(const char* fn);
        //! Open a bzip2-compressed file.
        void openBZ2(const std::string& fn)
        {
            openBZ2(fn.c_str());
        }
        //! Open a bzip2-compressed file.
        void openBZ2(const char* fn);

        void close()
        {
            reset();
        }
    private:
        typedef boost::iostreams::filtering_istream base_type;
    };

    class CompressedOutputStream:
        public boost::iostreams::filtering_ostream
    {
    public:
        static const int DEFAULT_COMPRESSION = 4;

        CompressedOutputStream() {}
        CompressedOutputStream(const char* fn, int clev = DEFAULT_COMPRESSION)
        {
            open(fn, clev);
        }
        CompressedOutputStream(const std::string& fn, int clev = DEFAULT_COMPRESSION)
        {
            open(fn, clev);
        }

        ~CompressedOutputStream();

        //! Open a gzip-compressed file.
        void open(const std::string& fn, int clev = DEFAULT_COMPRESSION)
        {
            open(fn.c_str(), clev);
        }
        //! Open a gzip-compressed file.
        void open(const char* fn, int clev = DEFAULT_COMPRESSION);
        //! Open a bzip2-compressed file.
        void openBZ2(const std::string& fn)
        {
            openBZ2(fn.c_str());
        }
        //! Open a bzip2-compressed file.
        void openBZ2(const char* fn);

        void close()
        {
            reset();
        }
    private:
        typedef boost::iostreams::filtering_ostream base_type;
    };

    //! Writes a '\\0' for false or '\\001' for true.
    void writeBinaryBool(std::ostream& out, bool val);

    //! Writes the integer out in msb -> lsb order.
    template <class IntType>
    void writeBinaryInt(std::ostream& out, IntType val)
    {
        for(int ii=sizeof(val)*8-8; ii>=0; ii-=8)
        {
            char ch = (val >> ii) & 0xff;
            out.put(ch);
        }
    }

    //! Writes the integer out in "zero-compressed" format, suitable for
    //! reading in by cgatools::util::readBinaryUIntZC(). See
    //! cgatools::util::readBinaryUIntZC().
    template <class IntType>
    void writeBinaryUIntZC(std::ostream& out, IntType val)
    {
        char bytes[sizeof(val)*2];
        char* ptr = bytes + sizeof(bytes);
        char flag = 0;
        for(;;)
        {
            char ch = (val & 0x7f) | flag;
            val >>= 7;
            if (0 == val)
            {
                *--ptr = ch;
                break;
            }
            flag = char(0x80);
            *--ptr = ch;
        }

        out.write(ptr, bytes+sizeof(bytes)-ptr);
    }

    //! Writes the string out in a binary format, suitable for reading
    //! in by cgatools::util::readBinaryString().
    void writeBinaryString(std::ostream& out, const std::string& val);

    //! Reads one byte interpreted as a boolean, such that 0 is false
    //! and anything else is true.
    void readBinaryBool(std::istream& in, bool* val);

    //! Reads the integer in msb -> lsb order.
    template <class IntType>
    void readBinaryInt(std::istream& in, IntType* pVal)
    {
        IntType& val = *pVal;
        val = 0;
        for(size_t ii=0; ii<sizeof(val); ii++)
        {
            char ch;
            in.get(ch);
            if (!in.good())
                throw Exception("failed to read binary int: unexpected eof");

            val <<= 8;
            val |= static_cast<uint8_t>(ch);
        }
    }
    
    //! Reads an integer stored in "zero-compressed" format. At each
    //! iteration, a byte is read in, and the lsb 7 bits are or'd into
    //! val. If the msb is set, then the val is shifted left by 7, and
    //! the next byte must be read in.
    template <class IntType>
    void readBinaryUIntZC(std::istream& in, IntType* pVal)
    {
        IntType& val = *pVal;
        IntType MAX_VAL = std::numeric_limits<IntType>::max() >> 7;
        val = 0;
        for(;;)
        {
            int ch = in.get();
            if (!in.good())
                throw Exception("failed to read zero-compressed binary int: unexpected eof");
            val |= ch & 0x7f;
            if ( 0 == (ch & 0x80) )
                break;
            if (val > MAX_VAL)
                throw Exception("failed to read zero-compressed binary int: overflow");
            val <<= 7;
        }
    }

    //! Reads the string written out by cgatools::util::writeBinaryString().
    void readBinaryString(std::istream& in, std::string* val);

} } // cgatools::util

#endif // CGATOOLS_UTIL_STREAMS_HPP_
