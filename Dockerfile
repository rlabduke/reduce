FROM alpine:latest

ARG REDUCE_VERSION

WORKDIR /app
RUN apk add --no-cache git cmake build-base && \
    git clone https://github.com/rlabduke/reduce.git src/reduce

WORKDIR /app/src/reduce
RUN git checkout ${REDUCE_VERSION}

WORKDIR /app/build/reduce
RUN cmake /app/src/reduce && make && make install

WORKDIR /app/data

ENTRYPOINT ["reduce"]
